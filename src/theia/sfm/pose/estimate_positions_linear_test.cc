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

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/pose/estimate_positions_linear.h"
#include "theia/sfm/pose/pairwise_translation_error.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/transformation/align_point_clouds.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/triplet_extractor.h"
#include "theia/util/map_util.h"
#include "theia/util/stringprintf.h"

namespace theia {

using Eigen::Vector3d;

namespace {

Camera RandomCamera() {
  Camera camera;
  camera.SetPosition(10 * Vector3d::Random());
  camera.SetOrientationFromAngleAxis(0.2 * Vector3d::Random());
  camera.SetImageSize(1000, 1000);
  camera.SetFocalLength(800);
  camera.SetAspectRatio(1.0);
  camera.SetSkew(0.0);
  camera.SetPrincipalPoint(500.0, 500.0);
  return camera;
}

Vector3d RelativeRotationFromTwoRotations(const Vector3d& rotation1,
                                          const Vector3d& rotation2,
                                          const double noise) {
  const Eigen::Matrix3d noisy_rotation =
      Eigen::AngleAxisd(DegToRad(noise), Vector3d::Random().normalized())
          .toRotationMatrix();

  Eigen::Matrix3d rotation_matrix1, rotation_matrix2;
  ceres::AngleAxisToRotationMatrix(rotation1.data(), rotation_matrix1.data());
  ceres::AngleAxisToRotationMatrix(rotation2.data(), rotation_matrix2.data());

  const Eigen::AngleAxisd relative_rotation(
      noisy_rotation * rotation_matrix2 * rotation_matrix1.transpose());
  return relative_rotation.angle() * relative_rotation.axis();
}

Vector3d RelativeTranslationFromTwoPositions(const Vector3d& position1,
                                             const Vector3d& position2,
                                             const Vector3d& rotation1,
                                             const double noise) {
  const Eigen::AngleAxisd noisy_translation(DegToRad(noise),
                                            Vector3d::Random().normalized());
  Eigen::Matrix3d rotation_matrix1;
  ceres::AngleAxisToRotationMatrix(rotation1.data(), rotation_matrix1.data());
  const Vector3d relative_translation =
        rotation_matrix1 * (position2 - position1).normalized();
  return noisy_translation * relative_translation;
}

// Aligns positions to the ground truth positions via a similarity
// transformation.
void AlignPositions(const std::unordered_map<ViewId, Vector3d>& gt_positions,
                    std::unordered_map<ViewId, Vector3d>* positions) {
  // Collect all positions into a vector.
  std::vector<Vector3d> gt_pos, pos;
  for (const auto& gt_position : gt_positions) {
    gt_pos.push_back(gt_position.second);
    const Vector3d& position = FindOrDie(*positions, gt_position.first);
    pos.push_back(position);
  }

  Eigen::Matrix3d rotation;
  Vector3d translation;
  double scale;
  AlignPointCloudsUmeyama(pos, gt_pos, &rotation, &translation, &scale);

  // Apply the similarity transformation.
  for (auto& position : *positions) {
    position.second = scale * (rotation * position.second) + translation;
  }
}

}  // namespace

class EstimatePositionsLinearTest : public ::testing::Test {
 public:
  void TestLinearPositionEstimator(const int num_views,
                                   const int num_tracks,
                                   const int num_view_pairs,
                                   const double pose_noise,
                                   const double position_tolerance) {
    // Set up the camera.
    SetupReconstruction(num_views, num_tracks);
    GetTwoViewInfos(num_view_pairs, pose_noise);
    GetTriplets();

    // Estimate the positions.
    LinearPositionEstimator position_estimator(options_,
                                               reconstruction_,
                                               triplets_);

    std::unordered_map<ViewId, Vector3d> estimated_positions;
    EXPECT_TRUE(position_estimator.EstimatePositions(orientations_,
                                                     &estimated_positions));
    EXPECT_EQ(estimated_positions.size(), positions_.size());

    // Align the positions and measure the error.
    AlignPositions(positions_, &estimated_positions);
    for (const auto& position : positions_) {
      const Vector3d& estimated_position =
          FindOrDie(estimated_positions, position.first);
      const double position_error =
          (position.second - estimated_position).norm();
      EXPECT_LT(position_error, position_tolerance)
          << "\ng.t. position = " << position.second.transpose()
          << "\nestimated position = " << estimated_position.transpose();
    }
  }

 protected:
  void SetUp() {
    srand(1234);
  }

  void SetupReconstruction(const int num_views, const int num_tracks) {
    // Create random views.
    std::vector<ViewId> view_ids;
    for (int i = 0; i < num_views; i++) {
      const ViewId view_id = reconstruction_.AddView(StringPrintf("%d", i));
      view_ids.push_back(view_id);

      // Create a random pose.
      *reconstruction_.MutableView(view_id)->MutableCamera() = RandomCamera();
      orientations_[view_id] =
          reconstruction_.View(view_id)->Camera().GetOrientationAsAngleAxis();
      positions_[view_id] =
          reconstruction_.View(view_id)->Camera().GetPosition();
    }

    // Add random tracks.
    for (int i = 0; i < num_tracks; i++) {
      // Shuffle the view ids so that we can obtain tracks in random views.
      std::random_shuffle(view_ids.begin(), view_ids.end());

      // Create a track that is seen in several views.
      Eigen::Vector4d point = Eigen::Vector4d::Random();
      point[2] += 20.0;
      point[3] = 1.0;
      std::vector<std::pair<ViewId, Feature> > features;
      for (int j = 0; j < view_ids.size(); j++) {
        const View* view = reconstruction_.View(view_ids[j]);
        Eigen::Vector2d pixel;
        view->Camera().ProjectPoint(point, &pixel);
        features.emplace_back(view_ids[j], pixel);
      }
      Track* track =
          reconstruction_.MutableTrack(reconstruction_.AddTrack(features));
      *track->MutablePoint() = point;
    }
  }

  void GetTwoViewInfos(const int num_view_pairs, const double pose_noise) {
    // Create a single connected component.
    std::vector<ViewId> view_ids;
    view_ids.push_back(0);
    for (int i = 1; i < positions_.size(); i++) {
      const ViewIdPair view_id_pair(i - 1, i);
      view_pairs_[view_id_pair] = CreateTwoViewInfo(view_id_pair, pose_noise);
      view_ids.push_back(i);
    }

    while (view_pairs_.size() < num_view_pairs) {
      std::random_shuffle(view_ids.begin(), view_ids.end());
      const ViewIdPair view_id_pair =
          (view_ids[0] < view_ids[1]) ? ViewIdPair(view_ids[0], view_ids[1])
                                      : ViewIdPair(view_ids[1], view_ids[0]);
      if (ContainsKey(view_pairs_, view_id_pair)) {
        continue;
      }

      view_pairs_[view_id_pair] = CreateTwoViewInfo(view_id_pair, pose_noise);
    }
  }

  TwoViewInfo CreateTwoViewInfo(const ViewIdPair& view_id_pair,
                                const double pose_noise) {
    CHECK_LT(view_id_pair.first, view_id_pair.second);
    TwoViewInfo info;
    info.focal_length_1 = 800.0;
    info.focal_length_2 = 800.0;

    // These objects will add noise to the relative pose.
    const Eigen::Vector2d noise = pose_noise * Eigen::Vector2d::Random();

    // Determine the relative rotation and add noise.
    info.rotation_2 = RelativeRotationFromTwoRotations(
        FindOrDie(orientations_, view_id_pair.first),
        FindOrDie(orientations_, view_id_pair.second),
        noise(0));

    // Determine the relative position and add noise.
    info.position_2 = RelativeTranslationFromTwoPositions(
        FindOrDie(positions_, view_id_pair.first),
        FindOrDie(positions_, view_id_pair.second),
        FindOrDie(orientations_, view_id_pair.first),
        noise(1));

    return info;
  }

  // Get the triplets
  void GetTriplets() {
    TripletExtractor extractor;
    std::vector<std::vector<ViewTriplet> > triplets_vec;
    CHECK(extractor.ExtractTripletsFromViewPairs(view_pairs_, &triplets_vec));
    triplets_ = triplets_vec[0];
  }

  LinearPositionEstimator::Options options_;
  std::unordered_map<ViewId, Vector3d> positions_;
  std::unordered_map<ViewId, Vector3d> orientations_;
  std::unordered_map<ViewIdPair, TwoViewInfo> view_pairs_;
  std::vector<ViewTriplet> triplets_;
  Reconstruction reconstruction_;
};

TEST_F(EstimatePositionsLinearTest, SmallTestNoNoise) {
  static const double kTolerance = 1e-4;
  static const int kNumViews = 4;
  static const int kNumTracksPerView = 50;
  static const int kNumViewPairs = 6;
  TestLinearPositionEstimator(kNumViews,
                              kNumTracksPerView,
                              kNumViewPairs,
                              0.0,
                              kTolerance);
}

TEST_F(EstimatePositionsLinearTest, SmallTestWithNoise) {
  static const double kTolerance = 0.1;
  static const int kNumViews = 4;
  static const int kNumTracksPerView = 50;
  static const int kNumViewPairs = 6;
  static const double kPoseNoiseDegrees = 1.0;
  TestLinearPositionEstimator(kNumViews,
                              kNumTracksPerView,
                              kNumViewPairs,
                              kPoseNoiseDegrees,
                              kTolerance);
}

}  // namespace theia
