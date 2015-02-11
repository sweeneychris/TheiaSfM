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

#ifndef THEIA_SFM_ESTIMATE_TWOVIEW_INFO_H_
#define THEIA_SFM_ESTIMATE_TWOVIEW_INFO_H_

#include <vector>
#include "theia/sfm/create_and_initialize_ransac_variant.h"

namespace theia {

struct CameraIntrinsics;
struct FeatureCorrespondence;
struct TwoViewInfo;

// Options for estimating two view infos.
struct EstimateTwoViewInfoOptions {
  // Type of Ransac variant to use.
  RansacType ransac_type = RansacType::RANSAC;

  // Maximum sampson error in pixels for correspondences to be inliers.
  double max_sampson_error_pixels = 4.0;

  // Ransac parameters.
  double expected_ransac_confidence = 0.9999;
  int min_ransac_iterations = 10;
  int max_ransac_iterations = 1000;
  bool use_mle = true;
};

// Estimates two view info for the given view pair from the correspondences. The
// correspondences should be in pixel coordinates (the method will perform
// normalization w.r.t focal length as necessary).
//
// There are three cases for estimating two view infos:
//   1) Both views are calibrated. In this case, the essential matrix is
//      estimated then decomposed to compute the two view info.
//   2) Both views are uncalibrated. The fundamental matrix is estimated and
//      decomposed to compute the two view info. NOTE: The quality of the focal
//      length is not always great with fundamental matrix decomposition.
//   3) One view is calibrated and one view is uncalibrated. NOTE: This case is
//      currently unsupported, and case 2) will be used instead.
//
// Returns true if a two view info could be successfully estimated and false if
// not.
bool EstimateTwoViewInfo(
    const EstimateTwoViewInfoOptions& options,
    const CameraIntrinsics& intrinsics1,
    const CameraIntrinsics& intrinsics2,
    const std::vector<FeatureCorrespondence>& correspondences,
    TwoViewInfo* twoview_info,
    std::vector<int>* inlier_indices);

}  // namespace theia

#endif  // THEIA_SFM_ESTIMATE_TWOVIEW_INFO_H_
