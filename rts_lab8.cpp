#include <iostream>
#include <atomic>
#include <vector>
#include <thread>
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>
#include <algorithm>
#include <limits>


// Структура координат
struct Position {
    double x;
    double y;
};

// Структура узла (Node)
template<typename T>
struct Node {
    typedef Node<T>* NodePtr;
    T data;
    NodePtr next;
    uint64_t version;

    Node(T val) : data(val), next(nullptr), version(0) {}
    Node() : next(nullptr), version(0) {}
};


template<typename T>
class LockFreeVersionedStack {
public:
    typedef Node<T>* NodePtr;
    typedef std::atomic<Node<T>*> AtomicNodePtr;
    typedef std::atomic<uint64_t> AtomicVersion;

    struct VersionedHead {
        AtomicVersion version;
        AtomicNodePtr head;
    };

    LockFreeVersionedStack(size_t readers_num) : subs_num_(readers_num) {
        stop_flag_.store(false);
        stack_.version.store(0);
        stack_.head.store(nullptr);

        subscribers_ = new std::atomic<uint64_t>[readers_num];
        for (size_t i = 0; i < readers_num; ++i) {
            subscribers_[i].store(0);
        }
    }

    ~LockFreeVersionedStack() {
        stop();
        NodePtr current = stack_.head.load();
        while (current) {
            NodePtr next = current->next;
            delete current;
            current = next;
        }
        for (auto node : trash_) {
            delete node;
        }
        delete[] subscribers_;
    }

    // Операция вставки 
    void push(T value) {
        NodePtr new_node = new Node<T>(value);
        new_node->next = stack_.head.load();
        new_node->version = stack_.version.load() + 1;
        stack_.head.store(new_node);
        stack_.version.fetch_add(1);
    }

    // Операция удаления 
    bool pop() {
        NodePtr old_node = stack_.head.load();
        if (old_node == nullptr) return false;

        NodePtr new_first_node = old_node->next;
        if (new_first_node != nullptr) {
            new_first_node->version = stack_.version.load() + 1;
        }

        stack_.head.store(new_first_node);
        stack_.version.fetch_add(1);

        update_trash(old_node);
        return true;
    }

    // Подписка читателя 
    bool subscribe(const unsigned int& id, NodePtr& stack_ptr, uint64_t last_seen_version) {
        // Сразу проверяем, не остановлена ли работа.
        // Если писатель остановился, новых версий не будет, выходим.
        if (stop_flag_.load()) return false;

        auto& sub = subscribers_[id];

        // Цикл ожидания новой версии.
        // Выполняется, если читатель уже на актуальной версии.
        while (stack_.version.load(std::memory_order_relaxed) <= last_seen_version) {
            if (stop_flag_.load()) return false;
        }

        uint64_t readed_version = sub.load();
        uint64_t current_version = stack_.version.load();

        // Попытка атомарно обновить версию подписчика
        while (!sub.compare_exchange_strong(readed_version, current_version)) {
            if (stop_flag_.load()) return false;
            current_version = stack_.version.load();
            readed_version = sub.load();
        }

        stack_ptr = stack_.head.load();
        return true;
    }

    void unsubscribe(const unsigned int& id) {
        subscribers_[id].store(0);
    }

    void stop() {
        stop_flag_.store(true);
    }

    bool is_stopped() {
        return stop_flag_.load();
    }

    uint64_t get_current_version() {
        return stack_.version.load();
    }

private:
    // Сборщик мусора 
    void update_trash(NodePtr old_node) {
        trash_.push_back(old_node);

        uint64_t min_version = std::numeric_limits<uint64_t>::max();
        bool has_readers = false;

        for (size_t i = 0; i < subs_num_; ++i) {
            uint64_t version = subscribers_[i].load();
            if (version == 0) continue;
            min_version = std::min(min_version, version);
            has_readers = true;
        }

        if (!has_readers) return;

        auto it = trash_.begin();
        while (it != trash_.end()) {
            if ((*it)->version < min_version) {
                delete* it;
                it = trash_.erase(it);
            }
            else {
                ++it;
            }
        }
    }

    VersionedHead stack_;
    std::atomic<uint64_t>* subscribers_;
    size_t subs_num_;
    std::vector<NodePtr> trash_;
    std::atomic<bool> stop_flag_;
};


double calculate_y(double x) {
    return -(x * x) + 4.0 * x;
}

// Проверка консистентности
bool check_consistency(const std::vector<Position>& data) {
    if (data.empty()) return true;

    double epsilon = std::numeric_limits<double>::epsilon();

    for (const auto& p : data) {
        double expected_y = calculate_y(p.x);
        // Строгая проверка с машинным эпсилоном
        if (std::abs(p.y - expected_y) > epsilon) return false;
    }

    for (size_t i = 0; i < data.size() - 1; ++i) {
        // Проверяем только убывание X
        if (data[i].x <= data[i + 1].x) return false;
    }
    return true;
}

void writer_thread_func(LockFreeVersionedStack<Position>& stack, double step) {
    double current_x = 0.0;
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> action_dist(0, 10);
    int log_counter = 0;

    std::cout << "[Writer] Начало генерации координат..." << std::endl;

    // Генерация последовательности
    while (current_x <= 4.0 && !stack.is_stopped()) {
        Position p;
        p.x = current_x;
        p.y = calculate_y(current_x);
        stack.push(p);

        // Логирование прогресса (будем выводить каждое 20-ое сообщение)
        if (++log_counter % 20 == 0) {
            std::cout << "[Writer] Выстрел: x=" << p.x << ", y=" << p.y
                << " (v: " << stack.get_current_version() << ")" << std::endl;
        }

        current_x += step;

        // Стресс-тест: удаление случайного кол-ва элементов
        if (action_dist(gen) > 8) {
            int pop_count = std::uniform_int_distribution<>(1, 3)(gen);
            for (int k = 0; k < pop_count; ++k) stack.pop();
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(7));
    }

    std::cout << "[Writer] Конец генерации. Ждём читающие потоки 1 секунду..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "[Writer] Остановка стека..." << std::endl;
    stack.stop();
}

void reader_thread_func(LockFreeVersionedStack<Position>& stack, int id, double step, std::atomic<int>& versions_read) {
    LockFreeVersionedStack<Position>::NodePtr current_head = nullptr;
    uint64_t last_seen_version = 0;

    // Чтение версий
    while (stack.subscribe(id, current_head, last_seen_version)) {
        last_seen_version = stack.get_current_version();

        if (current_head == nullptr) {
            stack.unsubscribe(id);
            continue;
        }

        std::vector<Position> snapshot;
        auto node = current_head;
        while (node != nullptr) {
            snapshot.push_back(node->data);
            node = node->next;
        }

        if (!check_consistency(snapshot)) {
            std::cerr << "[Reader " << id << "] ОШИБКА: Нарушена консистентность" << std::endl;
        }

        versions_read++;
        stack.unsubscribe(id);

        // Имитация тяжёлой работы, иначе они слишком быстро читать будут
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
}

int main() {
    setlocale(LC_ALL, "");
    const int NUM_READERS = 4;
    LockFreeVersionedStack<Position> stack(NUM_READERS);
    double step = pow(10, -3);

    std::cout << "Начало испытания \"По кому звонит Колокольчик\"..." << std::endl;

    std::vector<std::thread> readers;
    std::atomic<int> versions_read_counts[NUM_READERS];
    for (int i = 0; i < NUM_READERS; ++i) versions_read_counts[i] = 0;

    std::thread writer(writer_thread_func, std::ref(stack), step);

    for (int i = 0; i < NUM_READERS; ++i) {
        readers.emplace_back(reader_thread_func, std::ref(stack), i, step, std::ref(versions_read_counts[i]));
    }

    writer.join();
    for (auto& t : readers) {
        t.join();
    }

    std::cout << "\n--- ОТЧЁТ ---" << std::endl;
    std::cout << "Всего версий создано (Writer): " << stack.get_current_version() << std::endl;
    for (int i = 0; i < NUM_READERS; ++i) {
        std::cout << "Читатель " << i << " прочитал версий: " << versions_read_counts[i] << std::endl;
    }
}