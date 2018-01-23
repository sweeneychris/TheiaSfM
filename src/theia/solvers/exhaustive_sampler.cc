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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#include "theia/solvers/exhaustive_sampler.h"

#include <algorithm>
#include <glog/logging.h>
#include <vector>

#include "theia/solvers/sampler.h"

namespace theia {

ExhaustiveSampler::ExhaustiveSampler(
    const std::shared_ptr<RandomNumberGenerator>& rng,
    const int min_num_samples)
    : Sampler(rng, min_num_samples), i_(0), j_(1) {
  CHECK_EQ(this->min_num_samples_, 2) << "ExhaustiveSampler makes a hard "
                                         "assumption that the number of "
                                         "samples needed is 2.";
}

bool ExhaustiveSampler::Initialize(const int num_datapoints) {
  CHECK_GE(num_datapoints, this->min_num_samples_);
  num_datapoints_ = num_datapoints;
  return true;
}

// The next sample is determined deterministically.
bool ExhaustiveSampler::Sample(std::vector<int>* subset) {
  subset->emplace_back(i_);
  subset->emplace_back(j_);

  // Increment j and adjust the implicit for loops accordingly.
  ++j_;
  if (j_ >= num_datapoints_) {
    ++i_;
    // If i >= num_datapoints then we have enumerated all possible combinations.
    // We simply reset the outer loop (i) to 0 so that the combinations are
    // rengenerated.
    if (i_ >= num_datapoints_ - 1) {
      i_ = 0;
    }
    j_ = i_ + 1;
  }
  return true;
}

}  // namespace theia
