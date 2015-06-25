// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#include <glog/logging.h>
#include <functional>
#include "gtest/gtest.h"

#include "theia/util/mutable_priority_queue.h"

namespace theia {

TEST(MutablePriorityQueue, Constructor) {
  mutable_priority_queue<int, int> mpq;
  EXPECT_TRUE(mpq.empty());
  EXPECT_EQ(mpq.size(), 0);
}

TEST(MutablePriorityQueue, Clear) {
  mutable_priority_queue<int, int> mpq;
  mpq.insert(1, 0);
  mpq.insert(2, 1);
  EXPECT_EQ(mpq.size(), 2);
  mpq.clear();
  EXPECT_TRUE(mpq.empty());
}

TEST(MutablePriorityQueue, Erase) {
  mutable_priority_queue<int, int> mpq;
  mpq.insert(1, 1);
  mpq.insert(2, 2);
  mpq.insert(3, 3);

  // Remove (2, 2).
  mpq.erase(2);
  EXPECT_EQ(mpq.size(), 2);
  EXPECT_FALSE(mpq.contains(2));
}

TEST(MutablePriorityQueue, Pop) {
  mutable_priority_queue<int, int> mpq;
  mpq.insert(1, 1);
  mpq.insert(2, 2);
  mpq.insert(3, 3);

  // (1, 1) should be at the top.
  EXPECT_EQ(mpq.top().second, 1);
  mpq.pop();
  EXPECT_EQ(mpq.size(), 2);
  EXPECT_EQ(mpq.top().second, 2);
}

TEST(MutablePriorityQueue, Update) {
  mutable_priority_queue<int, int> mpq;
  mpq.insert(1, 1);
  mpq.insert(2, 2);
  mpq.insert(3, 3);

  // Make (2, 2) ---> (2, 0) which should push it to the top.
  mpq.update(2, 0);
  EXPECT_EQ(mpq.size(), 3);
  EXPECT_EQ(mpq.top().first, 2);
  EXPECT_EQ(mpq.top().second, 0);
}

TEST(MutablePriorityQueue, Find) {
  mutable_priority_queue<int, int> mpq;
  mpq.insert(1, 1);
  mpq.insert(2, 2);
  mpq.insert(3, 3);
  EXPECT_TRUE(mpq.contains(2));
  EXPECT_EQ(mpq.find(2), 2);
}

TEST(MutablePriorityQueue, MaxQueue) {
  mutable_priority_queue<int, int, std::less<int> > mpq;
  mpq.insert(1, 1);
  mpq.insert(2, 2);
  mpq.insert(3, 3);
  EXPECT_EQ(mpq.top().first, 3);
  EXPECT_EQ(mpq.top().second, 3);
}

}  // namespace theia
