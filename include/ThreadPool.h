// From: https://github.com/progschj/ThreadPool
// Modified to implement singleton pattern

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
  /**
   * @brief Singleton instance pointer
   */
  static inline ThreadPool *threadPool;

  /**
   * @brief Get ThreadPool instance
   * @param numThreads Number of threads to use
   * @return Reference to ThreadPool instance
   */
  static inline ThreadPool &getInstance(size_t numThreads = 4);

  /**
   * @brief Check if thread pool is running
   * @return Running status
   */
  inline bool isRunning() const;

  /**
   * @brief Get task queue size for specific thread
   * @param thread_id Thread ID to query
   * @return Queue size
   */
  inline size_t getQueueSize(size_t thread_id) const;

  /**
   * @brief Maximum allowed queue size per thread
   */
  static constexpr size_t MAX_QUEUE_SIZE = 1000;

  inline ~ThreadPool();

  /**
   * @brief Submit task to thread pool
   * @tparam F Function type
   * @tparam Args Argument types
   * @param f Function to execute
   * @param args Function arguments
   * @return Future object for task result
   */
  template <class F, class... Args>
  inline auto enqueue(F &&f, Args &&...args) -> std::future<typename std::invoke_result<F, Args...>::type>;

  /**
   * @brief Get number of threads in pool
   * @return Thread count
   */
  inline size_t getNumThreads() const;

protected:
  /**
   * @brief Constructor
   * @param threads Number of threads to create
   */
  inline ThreadPool(size_t threads);

  /**
   * @brief Set CPU affinity for thread
   * @param thread_id Thread identifier
   * @param cpu_id CPU core ID
   */
  inline void setAffinity(int cpu_id);

  // Prevent copying
  inline ThreadPool(const ThreadPool &)            = delete;
  inline ThreadPool &operator=(const ThreadPool &) = delete;

  // Worker threads
  std::vector<std::thread> workers;

  // Thread-local task queues and synchronization
  std::vector<std::queue<std::function<void()>>>        thread_local_tasks;
  std::vector<std::unique_ptr<std::mutex>>              thread_local_mutexes;
  std::vector<std::unique_ptr<std::condition_variable>> thread_local_conditions;
  std::atomic<size_t>                                   next_thread{0};
  std::atomic<bool>                                     stop{false};
  size_t                                                numThreads;

  // 修改任务计数器声明
  std::vector<std::atomic<size_t>> task_counts;

  // Singleton implementation
  static inline std::atomic<ThreadPool *> instance_{nullptr};
  static inline std::mutex                mutex_;

  /**
   * @brief Attempt to steal task from other threads
   * @param thief_id ID of stealing thread
   * @param task Reference to store stolen task
   * @return True if task was stolen
   */
  inline bool stealTask(size_t thief_id, std::function<void()> &task);

  /**
   * @brief Main worker thread function
   * @param thread_id Thread identifier
   */
  inline void workerFunction(size_t thread_id);
};

#include "ThreadPool.inl"

#endif /* THREADPOOL_H */