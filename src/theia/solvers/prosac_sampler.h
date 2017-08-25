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

#ifndef THEIA_SOLVERS_PROSAC_SAMPLER_H_
#define THEIA_SOLVERS_PROSAC_SAMPLER_H_

#include <algorithm>
#include <glog/logging.h>
#include <random>
#include <stdlib.h>
#include <vector>

#include "theia/solvers/sampler.h"
#include "theia/util/random.h"

namespace theia {
// Prosac sampler used for PROSAC implemented according to "Matching with PROSAC
// - Progressive Sampling Consensus" by Chum and Matas.
class ProsacSampler : public Sampler {
 public:
  ProsacSampler(const std::shared_ptr<RandomNumberGenerator>& rng,
                const int min_num_samples);
  ~ProsacSampler() {}

  bool Initialize(const int num_datapoints);
  // Set the sample such that you are sampling the kth prosac sample (Eq. 6).
  void SetSampleNumber(int k);

  // Samples the input variable data and fills the vector subset with the prosac
  // samples.
  // NOTE: This assumes that data is in sorted order by quality where data[i] is
  // of higher quality than data[j] for all i < j.
  bool Sample(std::vector<int>* subset_indices);

 private:
  int num_datapoints_;
  // Number of iterations of PROSAC before it just acts like ransac.
  int ransac_convergence_iterations_;

  // The kth sample of prosac sampling.
  int kth_sample_number_;
};

}  // namespace theia

#endif  // THEIA_SOLVERS_PROSAC_SAMPLER_H_
