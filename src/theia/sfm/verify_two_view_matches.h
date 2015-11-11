// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_VERIFY_TWO_VIEW_MATCHES_H_
#define THEIA_SFM_VERIFY_TWO_VIEW_MATCHES_H_

#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/twoview_info.h"

namespace theia {

struct VerifyTwoViewMatchesOptions {
  // Parameters for estimating the two view geometry.
  EstimateTwoViewInfoOptions estimate_twoview_info_options;

  // Minimum number of inlier matches in order to return true.
  int min_num_inlier_matches = 30;

  // Bundle adjust the two view geometry using inliers.
  bool bundle_adjustment = true;

  // If performing bundle adjustment, the 3D points are only considered inliers
  // if the initial triangulation error is less than this. This value is in
  // pixels.
  double triangulation_sq_max_reprojection_error = 15.0;
  // If performing bundle adjustment, the 3D points are only considered inliers
  // if the reprojection error after bundle adjustment is less than this. This
  // value is in pixels.
  double final_sq_max_reprojection_error = 5.0;
};

bool VerifyTwoViewMatches(
    const VerifyTwoViewMatchesOptions& options,
    const CameraIntrinsicsPrior& intrinsics1,
    const CameraIntrinsicsPrior& intrinsics2,
    const std::vector<FeatureCorrespondence>& correspondences,
    TwoViewInfo* twoview_info,
    std::vector<int>* inlier_indices);

}  // namespace theia

#endif  // THEIA_SFM_VERIFY_TWO_VIEW_MATCHES_H_
