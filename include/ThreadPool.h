// From: https://github.com/progschj/ThreadPool
// Modified to be a singleton class

#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

class ThreadPool {
public:
  static inline ThreadPool *threadPool;

  static inline ThreadPool &getInstance(size_t numThreads = 4);

  inline bool isRunning() const;

  inline size_t getQueueSize(size_t thread_id) const;

  static constexpr size_t MAX_QUEUE_SIZE = 1000;

  inline ~ThreadPool();

  template <class F, class... Args>
  inline auto enqueue(F &&f, Args &&...args) -> std::future<typename std::invoke_result<F, Args...>::type>;

  inline size_t getNumThreads() const;

protected:
  inline ThreadPool(size_t);

  inline void setAffinity(std::thread::id, int);

  inline ThreadPool(const ThreadPool &) = delete;

  inline ThreadPool &operator=(const ThreadPool &) = delete;

  // need to keep track of threads so we can join them
  std::vector<std::thread> workers;

  std::vector<std::queue<std::function<void()>>>        thread_local_tasks;
  std::vector<std::unique_ptr<std::mutex>>              thread_local_mutexes;
  std::vector<std::unique_ptr<std::condition_variable>> thread_local_conditions;
  std::atomic<size_t>                                   next_thread{0};
  std::atomic<bool>                                     stop{false};

  size_t numThreads;

  static inline std::atomic<ThreadPool *> instance_{nullptr};
  static inline std::mutex                mutex_;

  inline bool stealTask(size_t thief_id, std::function<void()> &task);

  inline void workerFunction(size_t thread_id);
};

#include "ThreadPool.inl"

#endif /* THREADPOOL_H */