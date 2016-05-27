// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include <glog/logging.h>
#include <Eigen/Core>
#include <vector>

#include "gtest/gtest.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/matching/indexed_feature_match.h"
#include "theia/matching/guided_epipolar_matcher.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/pose/util.h"
#include "theia/util/hash.h"
#include "theia/util/random.h"

namespace theia {

void TestGuidedEpipolarMatcher(const int num_valid_matches,
                               const int num_invalid_matches,
                               const int num_provided_matches) {
  static const int kNumDescriptorDimensions = 128;
  InitRandomGenerator();

  // Set up two cameras, with camera 1 being at the coordinate system origin.
  static const double kFocalLength = 800.0;
  static const double kPrincipalPoint = 1000.0;
  Camera camera1, camera2;
  camera1.SetFocalLength(kFocalLength);
  camera1.SetPrincipalPoint(kPrincipalPoint, kPrincipalPoint);
  camera2.SetFocalLength(kFocalLength);
  camera2.SetPrincipalPoint(kPrincipalPoint, kPrincipalPoint);
  camera2.SetOrientationFromRotationMatrix(ProjectToRotationMatrix(
      Eigen::Matrix3d::Identity() + 0.1 * Eigen::Matrix3d::Random()));
  camera2.SetPosition(Eigen::Vector3d::Random());

  // Create 3d points and reproject them into both images to form
  // correspondences.
  KeypointsAndDescriptors features1, features2;
  for (int i = 0; i < num_valid_matches; i++) {
    Eigen::Vector4d point(RandDouble(-2.0, 2.0),
                          RandDouble(-2.0, 2.0),
                          RandDouble(5.0, 10.0),
                          1.0);
    Eigen::Vector2d point1, point2;
    CHECK_GT(camera1.ProjectPoint(point, &point1), 0);
    CHECK_GT(camera2.ProjectPoint(point, &point2), 0);

    // Set the keypoints to correspond to the same 3D point.
    features1.keypoints.emplace_back(point1.x(),
                                     point1.y(),
                                     Keypoint::OTHER);
    features2.keypoints.emplace_back(point2.x(),
                                     point2.y(),
                                     Keypoint::OTHER);

    // Make the descriptors the same for each feature so that they will be
    // guaranteed to match.
    Eigen::VectorXf descriptor(kNumDescriptorDimensions);
    descriptor.setRandom();
    descriptor.normalize();
    features1.descriptors.emplace_back(descriptor);
    features2.descriptors.emplace_back(descriptor);
  }

  // Add bogus features to the image that have no matches.
  for (int i = 0; i < num_invalid_matches; i++) {
    features1.keypoints.emplace_back(
        RandDouble(0, 1000.0), RandDouble(0, 1000.0), Keypoint::OTHER);
    features2.keypoints.emplace_back(
        RandDouble(0, 1000.0), RandDouble(0, 1000.0), Keypoint::OTHER);
    features1.descriptors.emplace_back(
        Eigen::VectorXf::Random(kNumDescriptorDimensions).normalized());
    features2.descriptors.emplace_back(
        Eigen::VectorXf::Random(kNumDescriptorDimensions).normalized());
  }

  // Add some pre-computed matches if applicable.
  std::vector<IndexedFeatureMatch> matches;
  for (int i = 0; i < num_provided_matches; i++) {
    IndexedFeatureMatch match;
    match.feature1_ind = i;
    match.feature2_ind = i;
    match.distance =
        (features1.descriptors[i] - features2.descriptors[i]).norm();
    matches.emplace_back(match);
  }

  // Run guided matching.
  GuidedEpipolarMatcher::Options options;
  GuidedEpipolarMatcher matcher(options, camera1, camera2, features1,
                                features2);

  // Ensure that the guided matching returns true.
  EXPECT_TRUE(matcher.GetMatches(&matches));

  // Ensure enough matches were found. Due to the randomness of point generation
  // we set this to be 90% of the known number of valid matches.
  EXPECT_GT(matches.size(), 0.9 * num_valid_matches);

  // Ensure that all matches are valid matches.
  for (const IndexedFeatureMatch match : matches) {
    if (match.feature1_ind < num_valid_matches) {
      EXPECT_EQ(match.feature1_ind, match.feature2_ind);
    }
  }
}

TEST(GuidedEpipolarMatcherTest, NoInputMatchesSmall) {
  //TestGuidedEpipolarMatcher(1000, 500, 0);
  TestGuidedEpipolarMatcher(10, 0, 0);
}

TEST(GuidedEpipolarMatcherTest, NoInputMatchesLarge) {
  TestGuidedEpipolarMatcher(1000, 5000, 0);
}

TEST(GuidedEpipolarMatcherTest, WithInputMatchesSmall) {
  TestGuidedEpipolarMatcher(1000, 5000, 500);
}

TEST(GuidedEpipolarMatcherTest, WithInputMatchesLarge) {
  TestGuidedEpipolarMatcher(2000, 5000, 1000);
}

}  // namespace theia
