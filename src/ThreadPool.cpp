//
// Created by vujevic on 4/24/17.
//

#include "ThreadPool.hpp"
#include "Runnable.hpp"

Ref<Semaphore> createSemaphore(int value) {
    return Ref<Semaphore>(new Semaphore(value));
}

Ref<ThreadPool> createThreadPool(int numberOfThreads) {
    return Ref<ThreadPool>(new ThreadPool(numberOfThreads));
}

Semaphore::Semaphore(int value): value_(value) {

}

void Semaphore::post() {
    std::unique_lock<std::mutex> lock(mutex_);
    ++value_;
    condition_.notify_one();
}

void Semaphore::wait() {
    std::unique_lock<std::mutex> lock(mutex_);
    condition_.wait(lock, [&](){ return value_; });
    --value_;
}


ThreadPool::ThreadPool(int numOfThreads):
        numberOfThreads_(numOfThreads), terminated_(false) {

//    queueSem_ = createSemaphore(1);
//    activeSem_ = createSemaphore(0);

    for (int i =0; i < numOfThreads; i++) {
        threads_.emplace_back(std::thread(Worker(*this, i)));
    }
}

ThreadPool::~ThreadPool() {

    terminated_ = true;
}

void ThreadPool::executeAllAndWait(std::vector<Ref<Runnable>>& runnables) {

    {
        for (const auto& it : runnables) {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            tasks_.emplace(it);
        }
    }
    condition_.notify_one();
    terminated_ = true;
    condition_.notify_all();

    for (auto& it : threads_) it.join();

}
void Worker::operator()() {

    while(true) {
        Ref<Runnable> task;
        {
            std::unique_lock<std::mutex> lock(pool_.queue_mutex_);
            while(!pool_.terminated_ && pool_.tasks_.empty()) pool_.condition_.wait(lock);

            if (pool_.terminated_ && pool_.tasks_.empty()) return;

            task = pool_.tasks_.front();
            pool_.tasks_.pop();
        }
        task->run();
    }
}
