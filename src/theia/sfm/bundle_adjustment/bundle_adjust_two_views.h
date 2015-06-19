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

#ifndef THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUST_TWO_VIEWS_H_
#define THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUST_TWO_VIEWS_H_

#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/camera/camera.h"

namespace theia {

struct TwoViewInfo;

// Configuration parameters for two-view bundle adjustment. The two bools
// control whether intrinsics are additionally optimized which is useful for
// uncalibrated camera sets.
struct TwoViewBundleAdjustmentOptions {
  BundleAdjustmentOptions ba_options;
  bool constant_camera1_intrinsics = true;
  bool constant_camera2_intrinsics = true;
};

// Triangulates all 3d points and performs standard bundle adjustment on the
// points and cameras. The first camera is held constant so only the second
// camera and points are optimized. The cameras must be initialized with camera
// intrinsics. The camera pose should be initialized based on epipolar geometry
// estimation.
BundleAdjustmentSummary BundleAdjustTwoViews(
    const TwoViewBundleAdjustmentOptions& options,
    const std::vector<FeatureCorrespondence>& correspondence,
    Camera* camera1,
    Camera* camera2);

// Performs bundle adjustment to find the optimal rotation and translation
// describing the two views. This is done without the need for 3D points, as
// described in "Exact Two-Image Structure from Motion" by John Oliensis (PAMI
// 2002).
//
// NOTE: The correspondences must be normalized by the focal length and
// principal point.
BundleAdjustmentSummary BundleAdjustTwoViewsAngular(
    const BundleAdjustmentOptions& options,
    const std::vector<FeatureCorrespondence>& correspondences,
    TwoViewInfo* info);

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUST_TWO_VIEWS_H_
