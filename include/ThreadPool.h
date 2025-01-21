// From: https://github.com/progschj/ThreadPool
// Modified to be a singleton class

#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <utility>

class ThreadPool
{
public:
    static ThreadPool &getInstance(size_t);

    ~ThreadPool();

    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args)
        -> std::future<typename std::invoke_result<F, Args...>::type>;

    size_t getNumThreads() const { return numThreads; }

protected:
    ThreadPool(size_t);

    void setAffinity(std::thread::id, int);

    ThreadPool(const ThreadPool &) = delete;

    ThreadPool &operator=(const ThreadPool &) = delete;

    // need to keep track of threads so we can join them
    std::vector<std::thread> workers;

    std::vector<std::queue<std::function<void()>>> thread_local_tasks;
    std::vector<std::unique_ptr<std::mutex>> thread_local_mutexes;
    std::vector<std::unique_ptr<std::condition_variable>> thread_local_conditions;
    std::atomic<size_t> next_thread{0};
    bool stop;

    size_t numThreads;
};

inline ThreadPool &ThreadPool::getInstance(size_t numThreads = 4)
{
    static ThreadPool instance(numThreads);
    return instance;
}

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
    : stop(false), numThreads(threads)
{
    thread_local_tasks.resize(threads);
    for (size_t i = 0; i < threads; ++i)
    {
        thread_local_mutexes.push_back(std::make_unique<std::mutex>());
        thread_local_conditions.push_back(std::make_unique<std::condition_variable>());
    }

    for (size_t i = 0; i < threads; ++i)
        workers.emplace_back(
            [this, i]
            {
                setAffinity(std::this_thread::get_id(), i);

                while (true)
                {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(*thread_local_mutexes[i]);
                        thread_local_conditions[i]->wait(lock,
                                                         [this, i]
                                                         { return stop || !thread_local_tasks[i].empty(); });
                        if (stop && thread_local_tasks[i].empty())
                            return;
                        task = std::move(thread_local_tasks[i].front());
                        thread_local_tasks[i].pop();
                    }

                    task();
                }
            });
}

inline void ThreadPool::setAffinity(std::thread::id thread_id, int cpu_id)
{
    pthread_t thread = pthread_self();
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpu_id, &cpuset);
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
}

// add new work item to the pool
template <class F, class... Args>
auto ThreadPool::enqueue(F &&f, Args &&...args)
    -> std::future<typename std::invoke_result<F, Args...>::type>
{
    using return_type = typename std::invoke_result<F, Args...>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...));

    std::future<return_type> res = task->get_future();

    size_t target_thread = next_thread++ % workers.size();

    {
        std::unique_lock<std::mutex> lock(*thread_local_mutexes[target_thread]);

        // don't allow enqueueing after stopping the pool
        if (stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        thread_local_tasks[target_thread].emplace([task]()
                                                  { (*task)(); });
    }
    thread_local_conditions[target_thread]->notify_one();
    return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
    {
        for (size_t i = 0; i < thread_local_mutexes.size(); ++i)
        {
            std::unique_lock<std::mutex> lock(*thread_local_mutexes[i]);
            stop = true;
        }
    }
    for (auto &condition : thread_local_conditions)
    {
        condition->notify_one();
    }
    for (std::thread &worker : workers)
    {
        worker.join();
    }
}

#endif /* THREADPOOL_H */