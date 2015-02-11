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

#include "theia/sfm/verify_two_view_matches.h"

#include <glog/logging.h>
#include <vector>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/bundle_adjustment/bundle_adjust_two_views.h"
#include "theia/sfm/camera/camera_intrinsics.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/twoview_info.h"

namespace theia {

namespace {

bool BundleAdjustRelativePose(
    const std::vector<FeatureCorrespondence>& inliers,
    const CameraIntrinsics& intrinsics1,
    const CameraIntrinsics& intrinsics2,
    TwoViewInfo* info) {
  BundleAdjustmentOptions options;
  options.verbose = false;
  options.linear_solver_type = ceres::DENSE_SCHUR;

  // Perform two view bundle adjustment. Alternatively, we can try to optimize
  // the angular error of features
  BundleAdjustmentSummary summary =
      BundleAdjustTwoViews(options, inliers, intrinsics1, intrinsics2, info);
  return summary.success;
}

}  // namespace

bool VerifyTwoViewMatches(
    const VerifyTwoViewMatchesOptions& options,
    const CameraIntrinsics& intrinsics1,
    const CameraIntrinsics& intrinsics2,
    const std::vector<FeatureCorrespondence>& correspondences,
    TwoViewInfo* twoview_info,
    std::vector<int>* inlier_indices) {
  if (correspondences.size() < options.min_num_inlier_matches) {
    return false;
  }

  // Estimate the two view info. If we fail to estimate a two view info then do
  // not add this view pair to the verified matches.
  if (!EstimateTwoViewInfo(options.estimate_twoview_info_options,
                           intrinsics1,
                           intrinsics2,
                           correspondences,
                           twoview_info,
                           inlier_indices)) {
    return false;
  }

  // If there were not enough inliers, return false and do not bother to
  // (potentially) run bundle adjustment.
  if (inlier_indices->size() < options.min_num_inlier_matches) {
    return false;
  }

  // Bundle adjustment (optional).
  if (options.bundle_adjustment) {
    std::vector<FeatureCorrespondence> inliers(inlier_indices->size());
    for (int i = 0; i < inliers.size(); i++) {
      inliers[i] = correspondences[inlier_indices->at(i)];
    }
    if (!BundleAdjustRelativePose(inliers,
                                  intrinsics1,
                                  intrinsics2,
                                  twoview_info)) {
      return false;
    }
  }
  return true;
}

}  // namespace theia
