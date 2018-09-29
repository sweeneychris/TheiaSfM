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

#ifndef THEIA_SOLVERS_RANDOM_SAMPLER_H_
#define THEIA_SOLVERS_RANDOM_SAMPLER_H_

#include <algorithm>
#include <memory>
#include <numeric>
#include <stdlib.h>
#include <vector>

#include "theia/solvers/sampler.h"
#include "theia/util/random.h"

namespace theia {

// Random sampler used for RANSAC. This is guaranteed to generate a unique
// sample by performing a Fisher-Yates sampling.
class RandomSampler : public Sampler {
 public:
  RandomSampler(const std::shared_ptr<RandomNumberGenerator>& rng,
                const int min_num_samples);
  ~RandomSampler() {}

  bool Initialize(const int num_datapoints) override;

  // Samples the input variable data and fills the vector subset with the
  // random samples.
  bool Sample(std::vector<int>* subset_indices) override;

 private:
  std::vector<int> sample_indices_;
};

}  // namespace theia

#endif  // THEIA_SOLVERS_RANDOM_SAMPLER_H_
