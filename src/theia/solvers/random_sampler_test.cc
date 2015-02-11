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

#include <glog/logging.h>
#include <algorithm>
#include <vector>

#include "gtest/gtest.h"
#include "theia/solvers/random_sampler.h"

namespace theia {

namespace {

bool IsUnique(const std::vector<int>& vec) {
  std::vector<int> sorted_vec = vec;
  std::sort(sorted_vec.begin(), sorted_vec.end());
  return std::unique(sorted_vec.begin(), sorted_vec.end()) == sorted_vec.end();
}

}  // namespace

TEST(RandomSampler, UniqueMinimalSample) {
  static const int kMinNumSamples = 3;
  const std::vector<int> data_points = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  RandomSampler<int> sampler(kMinNumSamples);
  CHECK(sampler.Initialize());
  for (int i = 0; i < 100; i++) {
    std::vector<int> subset;
    EXPECT_TRUE(sampler.Sample(data_points, &subset));

    // Make sure that the sampling is unique.
    EXPECT_EQ(subset.size(), kMinNumSamples);
    EXPECT_TRUE(IsUnique(subset));
  }
}

}  // namespace theia
