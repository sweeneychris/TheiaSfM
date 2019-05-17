// Copyright (C) 2019 The Regents of the University of California (Regents).
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

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// May 2019

#ifndef THEIA_SFM_ESTIMATORS_ESTIMATE_RADIAL_DISTORTION_HOMOGRAPHY_H_
#define THEIA_SFM_ESTIMATORS_ESTIMATE_RADIAL_DISTORTION_HOMOGRAPHY_H_

#include <Eigen/Core>
#include <vector>

#include "theia/sfm/camera/division_undistortion_camera_model.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/pose/six_point_radial_distortion_homography.h"

namespace theia {

struct RansacParameters;
struct RansacSummary;

// imitates std::pair
// this is basically for radial fundamental matrix estimation
struct RadialDistortionFeatureCorrespondence {
 public:
  Feature feature_left;
  Feature feature_right;

  Feature normalized_feature_left;
  Feature normalized_feature_right;

  double focal_length_estimate_left = 1000.0;
  double focal_length_estimate_right = 1000.0;

  double min_radial_distortion = -5.0;
  double max_radial_distortion = 0.0;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

// Estimates the homography matrix from feature correspondences
// using the 6-pt radial distortion homography algorithm.
// Apart from the homography it also returns a distortion estimation for both
// cameras
bool EstimateRadialHomographyMatrix(
    const RansacParameters& ransac_params, const RansacType& ransac_type,
    const std::vector<RadialDistortionFeatureCorrespondence>&
        normalized_correspondences,
    RadialHomographyResult* result, RansacSummary* ransac_summary);

}  // namespace theia

#endif
