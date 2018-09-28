// Copyright (C) 2018 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#include <glog/logging.h>
#include <gtest/gtest.h>
#include <set>

#include "theia/math/reservoir_sampler.h"
#include "theia/util/random.h"

namespace theia {

TEST(ReservoirSampler, Sanity) {
  constexpr int kNumFeatures = 1000000;
  constexpr int kNumSampledFeatures = 1000;
  ReservoirSampler<int> reservoir_sampler(kNumSampledFeatures);

  // Add features ranging from 0 to 100. This should yield a random sample which
  // roughly corresponds to 0,1,....,99.
  for (int i = 0; i < kNumFeatures; i++) {
    reservoir_sampler.AddElementToSampler(i % kNumSampledFeatures);
  }

  RandomNumberGenerator rng;
  std::set<int> rand_samples;
  for (int i = 0; i < kNumSampledFeatures; i++) {
    rand_samples.insert(rng.RandInt(0, kNumSampledFeatures));
  }

  auto samples = reservoir_sampler.GetAllSamples();
  std::set<int> sorted_samples(samples.begin(), samples.end());
  LOG(INFO) << "Num unique reservoir samples: " << sorted_samples.size();
  LOG(INFO) << "Num unique random samples: " << rand_samples.size();

}

}  // namespace theia
