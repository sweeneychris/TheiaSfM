// Copyright (C) 2018 The Regents of the University of California (Regents).
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
// Author: Victor Fragoso (victor.fragoso@mail.wvu.edu)

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/estimators/camera_and_feature_correspondence_2d_3d.h"
#include "theia/sfm/estimators/estimate_rigid_transformation_2d_3d.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/sfm/rigid_transformation.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"
#include "theia/util/timer.h"

namespace theia {
namespace {
using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

static const int kNumCameras = 10;
static const int kNumPoints = 100;
static const double kFocalLength = 1000.0;
static const double kReprojectionError = 4.0;
static const double kAngularErrorThresh = 1.0;  // 1 degree.
static const int kSeed = 66;

inline Camera RandomCamera(RandomNumberGenerator* rng) {
  Camera camera;
  camera.SetPosition(rng->RandVector3d(-10.0, 10.0));
  camera.SetOrientationFromRotationMatrix(RandomRotation(10.0, rng));
  camera.SetImageSize(1000, 1000);
  camera.SetFocalLength(kFocalLength);
  camera.SetPrincipalPoint(kFocalLength / 2.0, kFocalLength / 2.0);
  return camera;
}

void ExecuteRandomTest(const RansacParameters& options,
                       const RigidTransformation& rigid_transformation,
                       const double inlier_ratio,
                       const double noise,
                       const double tolerance,
                       const int num_cameras,
                       RandomNumberGenerator* rng) {
  // Generate cameras.
  std::vector<Camera> cameras(num_cameras);
  for (int i = 0; i < cameras.size(); ++i) {
    cameras[i] = RandomCamera(rng);
  }
  
  // Create feature correspondences (inliers and outliers) and add noise if
  // appropriate.
  std::vector<CameraAndFeatureCorrespondence2D3D> correspondences;
  correspondences.reserve(kNumPoints);
  for (int i = 0; i < kNumPoints; i++) {
    CameraAndFeatureCorrespondence2D3D correspondence;

    // Set up random camera.
    correspondence.camera = cameras[i % cameras.size()];

    // Set up random 3D point and reproject it into the image. Make sure the
    // point is in front of the camera.
    double depth = -1;
    do {
      // Create a 3D point randomly, and try to ensure that it is in front of
      // the camera.
      correspondence.point3d = Eigen::Vector4d(rng->RandDouble(-5, 5),
                                               rng->RandDouble(-5, 5),
                                               rng->RandDouble(10, 20),
                                               1.0);

      depth = correspondence.camera.ProjectPoint(correspondence.point3d,
                                                 &correspondence.observation);
    } while (depth < 0);
    correspondences.emplace_back(std::move(correspondence));
  }

  // Add noise to the image observations.
  if (noise) {
    for (int i = 0; i < kNumPoints; ++i) {
      correspondences[i].observation += noise * rng->RandVector2d();
    }
  }

  // Add outliers.
  for (int i = 0; i < kNumPoints; ++i) {
    if (i > inlier_ratio * kNumPoints) {
      correspondences[i].observation = kFocalLength * rng->RandVector2d();
    }
  }

  // Apply the rigid transformation to the 3d points.
  for (int i = 0; i < kNumPoints; ++i) {
    const Eigen::Vector3d old_point = correspondences[i].point3d.hnormalized();
    // Since the 3d point is in front of the camera, we need to calculate the
    // final 3d point to estimate the rigid transformation.
    const Eigen::Vector3d new_point =
        rigid_transformation.rotation.transpose() *
        (old_point - rigid_transformation.translation);
    correspondences[i].point3d = new_point.homogeneous();
  }

  // Estimate the rigid transform.
  RigidTransformation estimated_rigid_transformation;
  RansacSummary ransac_summary;
  Timer timer;
  timer.Reset();
  EXPECT_TRUE(EstimateRigidTransformation2D3D(options,
                                              RansacType::RANSAC,
                                              correspondences,
                                              &estimated_rigid_transformation,
                                              &ransac_summary));
  const double elapsed_time = timer.ElapsedTimeInSeconds();

  VLOG(3) << "Ransac summary: \n Number of inliers: "
          << ransac_summary.inliers.size()
          << "\n Num. input data points: "
          << ransac_summary.num_input_data_points
          << "\n Num. iterations: "
          << ransac_summary.num_iterations
          << "\n Confidence: " << ransac_summary.confidence
          << "\n Time [sec]: " << elapsed_time;
  VLOG(3) << "Expected rotation: \n" << rigid_transformation.rotation
          << "\n Estimated rotation: \n"
          << estimated_rigid_transformation.rotation
          << "\n Expected translation: "
          << rigid_transformation.translation.transpose()
          << "\n Estimated translation: "
          << estimated_rigid_transformation.translation.transpose();

  // Expect that the inlier ratio is close to the ground truth.
  EXPECT_GT(static_cast<double>(ransac_summary.inliers.size()), 3);

  // Expect rigid transforms are close.
  const Eigen::Quaterniond gt_rotation(rigid_transformation.rotation);
  const Eigen::Quaterniond estimated_rotation(
      estimated_rigid_transformation.rotation);
  EXPECT_NEAR(RadToDeg(gt_rotation.angularDistance(estimated_rotation)),
              0.0, kAngularErrorThresh);
  EXPECT_NEAR((estimated_rigid_transformation.translation -
               rigid_transformation.translation).norm(),
              0.0, tolerance);
}

}  // namespace

class EstimateRigidTransformation : public ::testing::Test {
 public:
  static void SetUpTestCase() {
    CHECK_GT(kNumPoints, kNumCameras);
    rng = new RandomNumberGenerator(kSeed);
  }

  static void TearDownTestCase() {
    delete rng;
  }

  // TODO(vfragoso): Generate the cameras and points once.
  static RandomNumberGenerator* rng;
};

RandomNumberGenerator* EstimateRigidTransformation::rng = nullptr;

TEST_F(EstimateRigidTransformation, AllInliersNoNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(kSeed);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 1.0;
  const double kNoise = 0.0;
  const double kPoseTolerance = 1e-4;

  RigidTransformation rigid_transformation;
  rigid_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                        .normalized()).toRotationMatrix();
  rigid_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  ExecuteRandomTest(options,
                    rigid_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance,
                    kNumCameras,
                    rng);
}


TEST_F(EstimateRigidTransformation, AllInliersWithNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(kSeed);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 1.0;
  const double kNoise = 1.0;
  const double kPoseTolerance = 1e-2;

  RigidTransformation rigid_transformation;
  rigid_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                        .normalized()).toRotationMatrix();
  rigid_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  ExecuteRandomTest(options,
                    rigid_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance,
                    kNumCameras,
                    rng);
}


TEST_F(EstimateRigidTransformation, OutliersNoNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(kSeed);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 0.8;
  const double kNoise = 0.0;
  const double kPoseTolerance = 1e-4;

  RigidTransformation rigid_transformation;
  rigid_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                        .normalized()).toRotationMatrix();
  rigid_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  ExecuteRandomTest(options,
                    rigid_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance,
                    kNumCameras,
                    rng);
}


TEST_F(EstimateRigidTransformation, OutliersWithNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(kSeed);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 0.8;
  const double kNoise = 1.0;
  const double kPoseTolerance = 0.1;

  RigidTransformation rigid_transformation;
    rigid_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                        .normalized()).toRotationMatrix();
  rigid_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  ExecuteRandomTest(options,
                    rigid_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance,
                    kNumCameras,
                    rng);
}

// TODO(vfragoso): Add unit tests to localize a single camera.

}  // namespace theia
