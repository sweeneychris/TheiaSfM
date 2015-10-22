// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <bitset>
#include <string>
#include "gtest/gtest.h"

#include "theia/util/random.h"
#include "theia/matching/distance.h"

namespace theia {

namespace {

int kNumTrials = 100;

// Zero distance.
TEST(L2Distance, ZeroDistance) {
  Eigen::VectorXf descriptor1(4);
  descriptor1.setRandom();
  descriptor1.normalize();
  Eigen::VectorXf descriptor2 = descriptor1;
  L2 l2_dist;
  ASSERT_EQ(l2_dist(descriptor1, descriptor2), 0);
}

// Known distance.
TEST(L2Distance, KnownDistance) {
  InitRandomGenerator();
  const int num_dimensions = 128;
  Eigen::VectorXf descriptor1(num_dimensions);
  Eigen::VectorXf descriptor2(num_dimensions);
  for (int n = 0; n < kNumTrials; n++) {
    descriptor1.setRandom();
    descriptor2.setRandom();
    descriptor1.normalize();
    descriptor2.normalize();
    L2 l2_dist;
    const float dist =
        static_cast<float>(2.0 - 2.0 * descriptor1.dot(descriptor2));
    ASSERT_DOUBLE_EQ(l2_dist(descriptor1, descriptor2), dist);
  }
}

}  // namespace
}  // namespace theia
