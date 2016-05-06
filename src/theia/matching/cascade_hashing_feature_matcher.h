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

#ifndef THEIA_MATCHING_CASCADE_HASHING_FEATURE_MATCHER_H_
#define THEIA_MATCHING_CASCADE_HASHING_FEATURE_MATCHER_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "theia/matching/cascade_hasher.h"
#include "theia/matching/distance.h"
#include "theia/matching/feature_matcher.h"
#include "theia/util/hash.h"

namespace theia {

// Performs features matching between two sets of features using a cascade
// hashing approach. This hashing does not require any training and is extremely
// efficient but can only be used with float features like SIFT.
class CascadeHashingFeatureMatcher : public FeatureMatcher {
 public:
  explicit CascadeHashingFeatureMatcher(const FeatureMatcherOptions& options)
      : FeatureMatcher(options) {}
  ~CascadeHashingFeatureMatcher() {}

  // These methods are the same as the base class except that the HashedImage is
  // created as the descriptors are added.
  void AddImage(const std::string& image_name,
                const std::vector<Keypoint>& keypoints,
                const std::vector<Eigen::VectorXf>& descriptors) override;
  void AddImage(const std::string& image_name,
                const std::vector<Keypoint>& keypoints,
                const std::vector<Eigen::VectorXf>& descriptors,
                const CameraIntrinsicsPrior& intrinsics) override;
  void AddImage(const std::string& image_name) override;
  void AddImage(const std::string& image_name,
                const CameraIntrinsicsPrior& intrinsics) override;

 private:
  bool MatchImagePair(
      const KeypointsAndDescriptors& features1,
      const KeypointsAndDescriptors& features2,
      std::vector<FeatureCorrespondence>* matched_features) override;

  std::unordered_map<std::string, HashedImage> hashed_images_;
  std::unique_ptr<CascadeHasher> cascade_hasher_;

  DISALLOW_COPY_AND_ASSIGN(CascadeHashingFeatureMatcher);
};

}  // namespace theia

#endif  // THEIA_MATCHING_CASCADE_HASHING_FEATURE_MATCHER_H_
