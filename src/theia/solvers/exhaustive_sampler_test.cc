// Copyright (C) 2017 The Regents of the University of California (Regents).
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

#include <algorithm>
#include <glog/logging.h>
#include <memory>
#include <utility>
#include <vector>

#include "theia/solvers/exhaustive_sampler.h"
#include "theia/util/map_util.h"
#include "theia/util/random.h"
#include "gtest/gtest.h"

namespace theia {

TEST(ExhaustiveSampler, EnsureExhaustiveSample) {
  std::shared_ptr<RandomNumberGenerator> rng =
      std::make_shared<RandomNumberGenerator>(55);
  static const int kMinNumSamples = 2;
  static const int kNumDataPoints = 100;
  std::vector<int> data_points(kNumDataPoints);
  std::iota(data_points.begin(), data_points.end(), 0);

  ExhaustiveSampler sampler(rng, kMinNumSamples);
  CHECK(sampler.Initialize(data_points.size()));

  for (int i = 0; i < data_points.size(); i++) {
    for (int j = i + 1; j < data_points.size(); j++) {
      std::vector<int> subset;
      EXPECT_TRUE(sampler.Sample(&subset));

      // Make sure that the sampling is unique.
      EXPECT_EQ(subset.size(), kMinNumSamples);
      EXPECT_EQ(subset[0], i);
      EXPECT_EQ(subset[1], j);
    }
  }

  // The next sample after all combinations are enumerated should be (0, 1).
  std::vector<int> subset;
  EXPECT_TRUE(sampler.Sample(&subset));
  EXPECT_EQ(subset[0], 0);
  EXPECT_EQ(subset[1], 1);
}

}  // namespace theia
