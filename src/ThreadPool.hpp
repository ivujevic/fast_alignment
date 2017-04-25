//
// Created by vujevic on 4/24/17.
//

#pragma once
#include <iostream>

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>


#include "Ref.hpp"

class Runnable;
class Semaphore;
class ThreadPool;

class Worker {
public:
    Worker(ThreadPool& pool, int id):
            pool_(pool), id_(id) {};
    void operator()();

private:
    int id_;
    ThreadPool& pool_;
};


class Semaphore {
public:
    ~Semaphore() = default;
    int value() const { return value_; }

    void wait();
    void post();

    friend Ref<Semaphore> createSemaphore(int value);

private:

    Semaphore(int value);

    std::mutex mutex_;
    std::condition_variable condition_;
    int value_;
};

class ThreadPool {
public:
    ~ThreadPool();

    void executeAllAndWait(std::vector<Ref<Runnable>>& runnables);
    int getNumberOfThreads() const {
        return numberOfThreads_;
    }

    friend Ref<ThreadPool> createThreadPool(int num_threads);

private:

    ThreadPool(int numOfThreads);

    int numberOfThreads_;

    std::vector<std::thread> threads_;
    std::queue<Ref<Runnable>> tasks_;

    Ref<Semaphore> queueSem_;
    Ref<Semaphore> activeSem_;
    Ref<Semaphore> outSemaphore_;

    std::atomic<bool> terminated_;
    std::condition_variable condition_;
    std::mutex queue_mutex_;
    friend class Worker;
};

Ref<ThreadPool> createThreadPool(int a = 8);
