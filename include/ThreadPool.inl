#include "ThreadPool.h"

inline ThreadPool &ThreadPool::getInstance(size_t numThreads) {
  ThreadPool *instance = instance_.load(std::memory_order_acquire);
  if (instance == nullptr) {
    std::lock_guard<std::mutex> lock(mutex_);
    instance = instance_.load(std::memory_order_relaxed);
    if (instance == nullptr) {
      instance = new ThreadPool(numThreads);
      instance_.store(instance, std::memory_order_release);
    }
  }

  ThreadPool::threadPool = instance;

  return *instance;
}

inline bool ThreadPool::isRunning() const { return !stop.load(std::memory_order_acquire); }

inline size_t ThreadPool::getQueueSize(size_t thread_id) const {
  if (thread_id >= numThreads) return 0;
  std::lock_guard<std::mutex> lock(*thread_local_mutexes[thread_id]);
  return thread_local_tasks[thread_id].size();
}

inline size_t ThreadPool::getNumThreads() const { return numThreads; }

inline bool ThreadPool::stealTask(size_t thief_id, std::function<void()> &task) {
  for (size_t victim = 0; victim < numThreads; ++victim) {
    if (victim == thief_id) continue;
    std::lock_guard<std::mutex> lock(*thread_local_mutexes[victim]);
    if (!thread_local_tasks[victim].empty()) {
      task = std::move(thread_local_tasks[victim].front());
      thread_local_tasks[victim].pop();
      return true;
    }
  }
  return false;
}

inline void ThreadPool::workerFunction(size_t thread_id) {
  setAffinity(std::this_thread::get_id(), thread_id);

  while (true) {
    std::function<void()> task;
    bool                  has_task = false;

    {
      std::unique_lock<std::mutex> lock(*thread_local_mutexes[thread_id]);
      thread_local_conditions[thread_id]->wait(lock, [this, thread_id] {
        return stop || !thread_local_tasks[thread_id].empty();
      });

      if (!thread_local_tasks[thread_id].empty()) {
        task = std::move(thread_local_tasks[thread_id].front());
        thread_local_tasks[thread_id].pop();
        has_task = true;
      }
    }

    if (!has_task && !stop) { has_task = stealTask(thread_id, task); }

    if (!has_task && stop) { return; }

    if (has_task) { task(); }
  }
}

inline ThreadPool::ThreadPool(size_t threads) : stop(false), numThreads(threads) {
  thread_local_tasks.resize(threads);
  for (size_t i = 0; i < threads; ++i) {
    thread_local_mutexes.push_back(std::make_unique<std::mutex>());
    thread_local_conditions.push_back(std::make_unique<std::condition_variable>());
  }

  for (size_t i = 0; i < threads; ++i) workers.emplace_back([this, i] { workerFunction(i); });
}

inline void ThreadPool::setAffinity(std::thread::id thread_id, int cpu_id) {
  pthread_t thread = pthread_self();
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(cpu_id, &cpuset);
  pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
}

template <class F, class... Args>
auto ThreadPool::enqueue(F &&f, Args &&...args) -> std::future<typename std::invoke_result<F, Args...>::type> {
  if (stop.load(std::memory_order_acquire)) throw std::runtime_error("enqueue on stopped ThreadPool");

  using return_type = typename std::invoke_result<F, Args...>::type;

  auto task =
      std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

  std::future<return_type> res = task->get_future();

  size_t target_thread  = 0;
  size_t min_queue_size = getQueueSize(0);

  for (size_t i = 1; i < numThreads; ++i) {
    size_t current_size = getQueueSize(i);
    if (current_size < min_queue_size) {
      min_queue_size = current_size;
      target_thread  = i;
    }
  }

  if (min_queue_size >= MAX_QUEUE_SIZE) { throw std::runtime_error("Task queue is full"); }

  {
    std::unique_lock<std::mutex> lock(*thread_local_mutexes[target_thread]);
    thread_local_tasks[target_thread].emplace([task]() { (*task)(); });
  }

  thread_local_conditions[target_thread]->notify_one();
  return res;
}

inline ThreadPool::~ThreadPool() {
  stop.store(true, std::memory_order_release);

  for (auto &condition : thread_local_conditions) condition->notify_one();

  for (std::thread &worker : workers)
    if (worker.joinable()) worker.join();
}