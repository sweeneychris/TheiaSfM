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

#ifndef THEIA_SFM_MATCH_AND_VERIFY_FEATURES_H_
#define THEIA_SFM_MATCH_AND_VERIFY_FEATURES_H_

#include <Eigen/Core>
#include <vector>

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/matching/image_pair_match.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/verify_two_view_matches.h"

namespace theia {

struct CameraIntrinsics;

// The type of matching to perform.
enum class MatchingStrategy {
  BRUTE_FORCE = 0,
  CASCADE_HASHING = 1,
};

// Options for image matching with geometric verification.
struct MatchAndVerifyFeaturesOptions {
  // Number of threads for multithreading.
  int num_threads = 1;

  // Minimum number of inliers to consider the matches a good match.
  int min_num_inlier_matches = 30;

  // Matching strategy to use for establishing feature correspondences.
  MatchingStrategy matching_strategy;

  // Matching options for determining which feature matches are good matches.
  FeatureMatcherOptions feature_matcher_options;

  // Options for estimating the relative pose that is used for geometric
  // verification.
  VerifyTwoViewMatchesOptions geometric_verification_options;
};

// The two methods below perform image matching using the feature descriptors
// then perform geometric verification to determine the inlier
// correspondences. Geometric verification is done using the essential matrix or
// fundamental matrix depending on whether the image is calibrated. The
// CameraIntrinsics should be set to the default parameters
// (focal_length == 1.0)
// if the camera has not been explicitly calibrated.
bool MatchAndVerifyFeatures(
    const MatchAndVerifyFeaturesOptions& options,
    const std::vector<CameraIntrinsics>& intrinsics,
    const std::vector<std::vector<Keypoint> >& keypoints,
    const std::vector<std::vector<Eigen::VectorXf> >& descriptor,
    std::vector<ImagePairMatch>* matches);

bool MatchAndVerifyFeatures(
    const MatchAndVerifyFeaturesOptions& options,
    const std::vector<CameraIntrinsics>& intrinsics,
    const std::vector<std::vector<Keypoint> >& keypoints,
    const std::vector<std::vector<BinaryVectorX> >& descriptor,
    std::vector<ImagePairMatch>* matches);

}  // namespace theia

#endif  // THEIA_SFM_MATCH_AND_VERIFY_FEATURES_H_
