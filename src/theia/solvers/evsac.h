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
// Author: Victor Fragoso (vfragoso@cs.ucsb.edu)

#ifndef THEIA_SOLVERS_EVSAC_H_
#define THEIA_SOLVERS_EVSAC_H_

#include <statx/distributions/evd/common.h>

#include "theia/solvers/estimator.h"
#include "theia/solvers/evsac_sampler.h"
#include "theia/solvers/sample_consensus_estimator.h"

namespace theia {
// Estimate a model using EVSAC sampler.
template <class ModelEstimator>
class Evsac : public SampleConsensusEstimator<ModelEstimator> {
 public:
  typedef typename ModelEstimator::Datum Datum;
  typedef typename ModelEstimator::Model Model;

  // Params:
  // sorted_distances:  The matrix containing k L2 sorted distances. The
  //   matrix has num. of query features as rows and k columns.
  // predictor_threshold:  The threshold used to decide correct or
  //   incorrect matches/correspondences. The recommended value is 0.65.
  // fitting_method:  The fiting method MLE or QUANTILE_NLS (see statx doc).
  //   The recommended fitting method is the MLE estimation.
  Evsac(const RansacParameters& ransac_params,
        const ModelEstimator& estimator,
        const Eigen::MatrixXd& sorted_distances,
        const double predictor_threshold,
        const FittingMethod fitting_method)
      : SampleConsensusEstimator<ModelEstimator>(ransac_params, estimator),
        sorted_distances_(sorted_distances),
        predictor_threshold_(predictor_threshold),
        fitting_method_(fitting_method) {}

  ~Evsac() {}

  bool Initialize() {
    Sampler* prosac_sampler =
        new EvsacSampler<Datum>(this->estimator_.SampleSize(),
                                this->sorted_distances_,
                                this->predictor_threshold_,
                                this->fitting_method_);
    return SampleConsensusEstimator<ModelEstimator>::Initialize(prosac_sampler);
  }

 protected:
  // L2 descriptor sorted distances.
  // rows: num of reference features.
  // cols: k-th smallest distances.
  const Eigen::MatrixXd& sorted_distances_;
  // Threshold for predictor (MR-Rayleigh).
  const double predictor_threshold_;
  // The fitting method.
  const FittingMethod fitting_method_;

 private:
  DISALLOW_COPY_AND_ASSIGN(Evsac);
};
}  // namespace theia

#endif  // THEIA_SOLVERS_EVSAC_H_
