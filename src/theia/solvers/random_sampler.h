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

#include <stdlib.h>
#include <algorithm>
#include <vector>

#include "theia/solvers/sampler.h"
#include "theia/util/random.h"

namespace theia {

// Random sampler used for RANSAC. This is guaranteed to generate a unique
// sample by performing a Fisher-Yates sampling.
template <class Datum> class RandomSampler : public Sampler<Datum> {
 public:
  explicit RandomSampler(const int min_num_samples)
      : Sampler<Datum>(min_num_samples) {}
  ~RandomSampler() {}

  bool Initialize() {
    InitRandomGenerator();
    return true;
  }

  // Samples the input variable data and fills the vector subset with the
  // random samples.
  bool Sample(const std::vector<Datum>& data, std::vector<Datum>* subset) {
    subset->resize(this->min_num_samples_);
    std::vector<int> random_numbers(data.size());
    std::iota(random_numbers.begin(), random_numbers.end(), 0);

    for (int i = 0; i < this->min_num_samples_; i++) {
      std::swap(random_numbers[i], random_numbers[RandInt(i, data.size() - 1)]);
      (*subset)[i] = data[random_numbers[i]];
    }

    return true;
  }
};

}  // namespace theia

#endif  // THEIA_SOLVERS_RANDOM_SAMPLER_H_
