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

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/sfm/global_pose_estimation/robust_rotation_estimator.h"
#include "theia/sfm/transformation/align_rotations.h"
#include "theia/sfm/types.h"
#include "theia/util/map_util.h"
#include "theia/util/random.h"

namespace theia {

using Eigen::Vector3d;

namespace {

// Computes R_ij = R_j * R_i^t.
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

// return R_j = R_ij * R_i.
Vector3d ApplyRelativeRotation(const Vector3d& rotation1,
                               const Vector3d& relative_rotation) {
  Vector3d rotation2;
  Eigen::Matrix3d rotation1_matrix, relative_rotation_matrix;
  ceres::AngleAxisToRotationMatrix(
      rotation1.data(), ceres::ColumnMajorAdapter3x3(rotation1_matrix.data()));
  ceres::AngleAxisToRotationMatrix(
      relative_rotation.data(),
      ceres::ColumnMajorAdapter3x3(relative_rotation_matrix.data()));

  const Eigen::Matrix3d rotation2_matrix =
      relative_rotation_matrix * rotation1_matrix;
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(rotation2_matrix.data()), rotation2.data());
  return rotation2;
}

// Aligns rotations to the ground truth rotations via a similarity
// transformation.
void AlignOrientations(const std::unordered_map<ViewId, Vector3d>& gt_rotations,
                       std::unordered_map<ViewId, Vector3d>* rotations) {
  // Collect all rotations into a vector.
  std::vector<Vector3d> gt_rot, rot;
  std::unordered_map<int, int> index_to_view_id;
  int current_index = 0;
  for (const auto& gt_rotation : gt_rotations) {
    gt_rot.emplace_back(gt_rotation.second);
    rot.emplace_back(FindOrDie(*rotations, gt_rotation.first));

    index_to_view_id[current_index] = gt_rotation.first;
    ++current_index;
  }

  AlignRotations(gt_rot, &rot);

  for (int i = 0; i < rot.size(); i++) {
    const ViewId view_id = FindOrDie(index_to_view_id, i);
    (*rotations)[view_id] = rot[i];
  }
}

}  // namespace

class EstimateRotationsRobustTest : public ::testing::Test {
 public:
  void TestRobustRotationEstimator(const int num_views,
                                   const int num_view_pairs,
                                   const double rotation_noise,
                                   const double rotation_tolerance_degrees) {
    // Set up the camera.
    CreateGTOrientations(num_views);
    GetRelativeRotations(num_view_pairs, rotation_noise);

    // Estimate the rotations.
    RobustRotationEstimator::Options options;
    RobustRotationEstimator rotation_estimator(options);

    // Set the initial rotation estimations.
    std::unordered_map<ViewId, Vector3d> estimated_rotations;
    InitializeRotationsFromSpanningTree(&estimated_rotations);

    EXPECT_TRUE(rotation_estimator.EstimateRotations(view_pairs_,
                                                     &estimated_rotations));
    EXPECT_EQ(estimated_rotations.size(), orientations_.size());

    // Align the rotations and measure the error.
    AlignOrientations(orientations_, &estimated_rotations);
    for (const auto& rotation : orientations_) {
      const Vector3d& estimated_rotation =
          FindOrDie(estimated_rotations, rotation.first);
      const Vector3d relative_rotation = RelativeRotationFromTwoRotations(
          estimated_rotation, rotation.second, 0.0);
      const double angular_error = RadToDeg(relative_rotation.norm());

      EXPECT_LT(angular_error, rotation_tolerance_degrees)
          << "\ng.t. rotations = " << rotation.second.transpose()
          << "\nestimated rotations = " << estimated_rotation.transpose();
    }
  }

 protected:
  void SetUp() {
    srand(4567);
  }

  void CreateGTOrientations(const int num_views) {
    static const double kRotationScale = 0.2;
    // Create random orientations.
    for (int i = 0; i < num_views; i++) {
      orientations_[i] = kRotationScale * Vector3d::Random();
    }
  }

  void GetRelativeRotations(const int num_view_pairs, const double pose_noise) {
    // Create a set of view id pairs that will contain a spanning tree.
    for (int i = 1; i < orientations_.size(); i++) {
      const ViewIdPair view_id_pair(i - 1, i);
      view_pairs_[view_id_pair].rotation_2 = RelativeRotationFromTwoRotations(
          FindOrDie(orientations_, view_id_pair.first),
          FindOrDie(orientations_, view_id_pair.second),
          pose_noise);
    }

    // Add random edges.
    InitRandomGenerator();
    while (view_pairs_.size() < num_view_pairs) {
      ViewIdPair view_id_pair(RandInt(0, orientations_.size() - 1),
                              RandInt(0, orientations_.size() - 1));
      // Ensure the first id is smaller than the second id.
      if (view_id_pair.first > view_id_pair.second) {
        view_id_pair = ViewIdPair(view_id_pair.second, view_id_pair.first);
      }

      // Do not add the view pair if it already exists.
      if (view_id_pair.first == view_id_pair.second ||
          ContainsKey(view_pairs_, view_id_pair)) {
        continue;
      }

      view_pairs_[view_id_pair].rotation_2 = RelativeRotationFromTwoRotations(
          FindOrDie(orientations_, view_id_pair.first),
          FindOrDie(orientations_, view_id_pair.second),
          pose_noise);
    }
  }

  // Initialize the rotations from a spanning tree.
  void InitializeRotationsFromSpanningTree(
      std::unordered_map<ViewId, Vector3d>* initial_orientations) {
    // Set the first view to be at the origin.
    (*initial_orientations)[0] = Vector3d::Zero();
    for (int i = 1; i < orientations_.size(); i++) {
      (*initial_orientations)[i] = ApplyRelativeRotation(
          FindOrDie(*initial_orientations, i - 1),
          FindOrDieNoPrint(view_pairs_, ViewIdPair(i - 1, i)).rotation_2);
    }
  }

  std::unordered_map<ViewId, Vector3d> orientations_;
  std::unordered_map<ViewIdPair, TwoViewInfo> view_pairs_;
};

TEST_F(EstimateRotationsRobustTest, SmallTestNoNoise) {
  static const double kTolerance = 1e-8;
  static const int kNumViews = 4;
  static const int kNumViewPairs = 6;
  TestRobustRotationEstimator(kNumViews, kNumViewPairs, 0.0, kTolerance);
}

TEST_F(EstimateRotationsRobustTest, SmallTestWithNoise) {
  static const double kToleranceDegrees = 1.0;
  static const int kNumViews = 4;
  static const int kNumViewPairs = 6;
  static const double kPoseNoiseDegrees = 1.0;
  TestRobustRotationEstimator(kNumViews,
                              kNumViewPairs,
                              kPoseNoiseDegrees,
                              kToleranceDegrees);
}

TEST_F(EstimateRotationsRobustTest, LargeTestWithNoise) {
  static const double kToleranceDegrees = 5.0;
  static const int kNumViews = 100;
  static const int kNumViewPairs = 800;
  static const double kPoseNoiseDegrees = 2.0;
  TestRobustRotationEstimator(kNumViews,
                              kNumViewPairs,
                              kPoseNoiseDegrees,
                              kToleranceDegrees);
}

}  // namespace theia
