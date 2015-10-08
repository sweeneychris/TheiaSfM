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

#ifndef THEIA_SFM_GLOBAL_POSE_ESTIMATION_NONLINEAR_ROTATION_ESTIMATOR_H_
#define THEIA_SFM_GLOBAL_POSE_ESTIMATION_NONLINEAR_ROTATION_ESTIMATOR_H_

#include <Eigen/Core>
#include <unordered_map>

#include "theia/sfm/global_pose_estimation/rotation_estimator.h"
#include "theia/sfm/types.h"
#include "theia/util/hash.h"

namespace theia {

// Computes the global rotations given relative rotations and an initial guess
// for the global orientations. Nonlinear optimization is performed with Ceres
// using a SoftL1 loss function to be robust to outliers.
class NonlinearRotationEstimator : public RotationEstimator {
 public:
  NonlinearRotationEstimator() : robust_loss_width_(0.1) {}
  explicit NonlinearRotationEstimator(const double robust_loss_width)
      : robust_loss_width_(robust_loss_width) {}

  // Estimates the global orientations of all views based on an initial
  // guess. Returns true on successful estimation and false otherwise.
  bool EstimateRotations(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations);

 private:
  const double robust_loss_width_;
};

}  // namespace theia

#endif  // THEIA_SFM_GLOBAL_POSE_ESTIMATION_NONLINEAR_ROTATION_ESTIMATOR_H_
