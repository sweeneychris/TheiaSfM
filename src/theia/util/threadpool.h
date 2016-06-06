// Copyright (C) 2014 The Regents of the University of California (Regents).
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

// This implementation is based on the implementation found here:
// https://github.com/progschj/ThreadPool
//
// As such, we include the license here.

// Copyright (c) 2012 Jakob Progsch, Václav Zeman
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
//    1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
//
//    2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
//
//    3. This notice may not be removed or altered from any source
//    distribution.

#ifndef THEIA_UTIL_THREADPOOL_H_
#define THEIA_UTIL_THREADPOOL_H_

#include <glog/logging.h>

#include <chrono>
#include <condition_variable>
#include <deque>
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

#include "theia/util/util.h"

namespace theia {

// A simple threadpool implementation.
class ThreadPool {
 public:
  // All the threads are created upon construction.
  explicit ThreadPool(const int num_threads);
  ~ThreadPool();

  // Adds a task to the threadpool.
  template <class F, class... Args>
  auto Add(F&& f, Args&& ... args)
      ->std::future<typename std::result_of<F(Args...)>::type>;

  // Waits for all current tasks/threads to finish.
  void WaitForTasksToFinish();

 private:
  // Keep track of threads so we can join them
  std::vector<std::thread> workers;
  // The task queue
  std::queue<std::function<void()> > tasks;

  // Synchronization
  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;

  DISALLOW_COPY_AND_ASSIGN(ThreadPool);
};

// add new work item to the pool
template <class F, class... Args>
auto ThreadPool::Add(F&& f, Args&& ... args)
    ->std::future<typename std::result_of<F(Args...)>::type> {
  using return_type = typename std::result_of<F(Args...)>::type;

  auto task = std::make_shared<std::packaged_task<return_type()> >(
      std::bind(std::forward<F>(f), std::forward<Args>(args)...));

  std::future<return_type> res = task->get_future();
  {
    std::unique_lock<std::mutex> lock(queue_mutex);

    // don't allow enqueueing after stopping the pool
    CHECK(!stop) << "The ThreadPool object has been destroyed! Cannot add more "
                    "tasks to the ThreadPool!";

    tasks.emplace([task]() {
      (*task)();
    });
  }
  condition.notify_one();
  return res;
}

}  // namespace theia

#endif  // THEIA_UTIL_THREADPOOL_H_
