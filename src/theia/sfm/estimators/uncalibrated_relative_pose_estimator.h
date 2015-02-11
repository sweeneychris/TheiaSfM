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

#ifndef THEIA_SFM_ESTIMATORS_UNCALIBRATED_RELATIVE_POSE_ESTIMATOR_H_
#define THEIA_SFM_ESTIMATORS_UNCALIBRATED_RELATIVE_POSE_ESTIMATOR_H_

#include <Eigen/Core>
#include <vector>

#include "theia/solvers/estimator.h"
#include "theia/util/util.h"
#include "theia/matching/feature_correspondence.h"

namespace theia {

// Relative pose information computed from two uncalibrated views. It is assumed
// that the first view has identity rotation and a position at the origin.
struct UncalibratedRelativePose {
  Eigen::Matrix3d fundamental_matrix;
  double focal_length1;
  double focal_length2;
  Eigen::Matrix3d rotation;
  Eigen::Vector3d position;
};

// An estimator for computing the relative pose from 8 feature correspondences
// (via decomposition of the fundamental matrix).
//
// NOTE: Feature correspondences must be in pixel coordinates with the principal
// point removed i.e. principal point at (0, 0). This also assumes negligible
// skew (which is reasonable for most cameras).
class UncalibratedRelativePoseEstimator
    : public Estimator<FeatureCorrespondence, UncalibratedRelativePose> {
 public:
  UncalibratedRelativePoseEstimator() {}

  // 8 correspondences are needed to determine a fundamental matrix and thus a
  // relative pose.
  double SampleSize() const { return 8; }

  // Estimates candidate relative poses from correspondences.
  bool EstimateModel(
      const std::vector<FeatureCorrespondence>& centered_correspondences,
      std::vector<UncalibratedRelativePose>* relative_poses) const;

  // The error for a correspondences given a model. This is the squared sampson
  // error.
  double Error(const FeatureCorrespondence& centered_correspondence,
               const UncalibratedRelativePose& relative_poses) const;

 private:
  DISALLOW_COPY_AND_ASSIGN(UncalibratedRelativePoseEstimator);
};

}  // namespace theia

#endif  // THEIA_SFM_ESTIMATORS_UNCALIBRATED_RELATIVE_POSE_ESTIMATOR_H_
