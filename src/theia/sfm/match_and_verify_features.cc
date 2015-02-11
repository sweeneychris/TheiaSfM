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

#include "theia/sfm/match_and_verify_features.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <mutex>
#include <vector>

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/matching/brute_force_feature_matcher.h"
#include "theia/matching/cascade_hashing_feature_matcher.h"
#include "theia/matching/distance.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/matching/feature_matcher.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/matching/image_pair_match.h"
#include "theia/sfm/camera/camera_intrinsics.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/verify_two_view_matches.h"
#include "theia/util/threadpool.h"

namespace theia {
namespace {

// Templated method so that we can perform the same operations of binary and
// float descriptors in one method.
template <class DistanceMetric, class DescriptorType>
bool MatchAndVerifyFeatures(
    const MatchAndVerifyFeaturesOptions& options,
    const std::vector<CameraIntrinsics>& intrinsics,
    const std::vector<std::vector<Keypoint>*>& keypoints,
    const std::vector<std::vector<DescriptorType>*>& descriptors,
    FeatureMatcher<DistanceMetric>* matcher,
    std::vector<ImagePairMatch>* matches) {
  CHECK_GT(options.num_threads, 0);
  CHECK_GT(options.min_num_inlier_matches, 0);
  CHECK_EQ(intrinsics.size(), keypoints.size());
  CHECK_EQ(keypoints.size(), descriptors.size());
  CHECK_NOTNULL(matches);

  // Match features.
  std::vector<ImagePairMatch> image_pair_matches;
  for (int i = 0; i < keypoints.size(); i++) {
    matcher->AddImage(keypoints[i], descriptors[i], intrinsics[i]);
  }

  // Set the options to ensure they are consistent with each other.
  FeatureMatcherOptions matcher_options = options.feature_matcher_options;
  matcher_options.num_threads = options.num_threads;
  matcher_options.min_num_feature_matches = options.min_num_inlier_matches;
  VerifyTwoViewMatchesOptions verification_options =
      options.geometric_verification_options;
  verification_options.min_num_inlier_matches = options.min_num_inlier_matches;

  matcher->MatchImagesWithGeometricVerification(matcher_options,
                                                verification_options,
                                                matches);
  return true;
}

}  // namespace

bool MatchAndVerifyFeatures(
    const MatchAndVerifyFeaturesOptions& options,
    const std::vector<CameraIntrinsics>& intrinsics,
    const std::vector<std::vector<Keypoint>*>& keypoints,
    const std::vector<std::vector<Eigen::VectorXf>*>& descriptor,
    std::vector<ImagePairMatch>* matches) {
  std::unique_ptr<FeatureMatcher<L2>> matcher;
  if (options.matching_strategy == MatchingStrategy::CASCADE_HASHING) {
    matcher.reset(new CascadeHashingFeatureMatcher);
  } else if (options.matching_strategy == MatchingStrategy::BRUTE_FORCE) {
    matcher.reset(new BruteForceFeatureMatcher<L2>);
  } else {
    LOG(FATAL) << "Invalid matching strategy specified.";
  }

  return MatchAndVerifyFeatures(options,
                                intrinsics,
                                keypoints,
                                descriptor,
                                matcher.get(),
                                matches);
}

bool MatchAndVerifyFeatures(
    const MatchAndVerifyFeaturesOptions& options,
    const std::vector<CameraIntrinsics>& intrinsics,
    const std::vector<std::vector<Keypoint>*>& keypoints,
    const std::vector<std::vector<BinaryVectorX>*>& descriptor,
    std::vector<ImagePairMatch>* matches) {
  std::unique_ptr<FeatureMatcher<Hamming>> matcher;
  if (options.matching_strategy == MatchingStrategy::CASCADE_HASHING) {
    LOG(WARNING) << "Invalid matching specified for binary vectors. Using "
        "brute force matching instead.";
    matcher.reset(new BruteForceFeatureMatcher<Hamming>);
  } else if (options.matching_strategy == MatchingStrategy::BRUTE_FORCE) {
    matcher.reset(new BruteForceFeatureMatcher<Hamming>);
  } else {
    LOG(FATAL) << "Invalid matching strategy specified.";
  }

  return MatchAndVerifyFeatures(options,
                                intrinsics,
                                keypoints,
                                descriptor,
                                matcher.get(),
                                matches);
}

}  // namespace theia
