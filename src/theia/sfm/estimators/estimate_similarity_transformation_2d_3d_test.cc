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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include "theia/math/util.h"
#include "theia/sfm/estimators/estimate_similarity_transformation_2d_3d.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/sfm/similarity_transformation.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"

namespace theia {

static const double kFocalLength = 1000.0;
static const double kReprojectionError = 4.0;
static const int kNumCameras = 100;

Camera RandomCamera() {
  Camera camera;
  camera.SetPosition(10 * Eigen::Vector3d::Random());
  camera.SetOrientationFromAngleAxis(0.2 * Eigen::Vector3d::Random());
  camera.SetImageSize(1000, 1000);
  camera.SetFocalLength(kFocalLength);
  camera.SetAspectRatio(1.0);
  camera.SetSkew(0.0);
  camera.SetPrincipalPoint(kFocalLength / 2.0, kFocalLength / 2.0);
  return camera;
}

void ExecuteRandomTest(
    const RansacParameters& params,
    const SimilarityTransformation& similarity_transformation,
    const double inlier_ratio,
    const double noise,
    const double tolerance) {
  std::vector<CameraAndFeatureCorrespondence2D3D> correspondences;
  for (int i = 0; i < kNumCameras; i++) {
    CameraAndFeatureCorrespondence2D3D correspondence;

    // Set up random camera.
    correspondence.camera = RandomCamera();

    // Set up random 3D point and reproject it into the image. Make sure the
    // point is in front of the camera.
    double depth = -1;
    do {
      // Create a 3D point randomly, and try to ensure that it is in front of
      // the camera.
      correspondence.point3d = Eigen::Vector4d(RandDouble(-5, 5),
                                               RandDouble(-5, 5),
                                               RandDouble(10, 20),
                                               1.0);

      depth = correspondence.camera.ProjectPoint(correspondence.point3d,
                                                 &correspondence.observation);
    } while (depth < 0);
    correspondences.emplace_back(correspondence);
  }

  // Add noise to the image observations.
  if (noise) {
    for (int i = 0; i < kNumCameras; i++) {
      correspondences[i].observation += noise * Eigen::Vector2d::Random();
    }
  }

  // Add outliers.
  for (int i = 0; i < kNumCameras; i++) {
    if (i > inlier_ratio * kNumCameras) {
      correspondences[i].observation = kFocalLength * Eigen::Vector2d::Random();
    }
  }

  // Apply the similarity transformation to the 3d points.
  for (int i = 0; i < kNumCameras; i++) {
    const Eigen::Vector3d old_point = correspondences[i].point3d.hnormalized();
    const Eigen::Vector3d new_point = similarity_transformation.scale *
                                          similarity_transformation.rotation *
                                          old_point +
                                      similarity_transformation.translation;
    correspondences[i].point3d = new_point.homogeneous();
  }

  // Estimate the similarity transformation;
  SimilarityTransformation estimated_similarity_transformation;
  RansacSummary ransac_summary;
  EXPECT_TRUE(EstimateSimilarityTransformation2D3D(
      params,
      RansacType::RANSAC,
      correspondences,
      &estimated_similarity_transformation,
      &ransac_summary));

  // We should have found at least one good solution.
  EXPECT_GT(ransac_summary.inliers.size(), 4);

  // Expect poses are near.
  EXPECT_TRUE(test::ArraysEqualUpToScale(
      9,
      estimated_similarity_transformation.rotation.data(),
      similarity_transformation.rotation.data(),
      tolerance));
  EXPECT_TRUE(test::ArraysEqualUpToScale(
      3,
      estimated_similarity_transformation.translation.data(),
      similarity_transformation.translation.data(),
      tolerance));
  EXPECT_NEAR(estimated_similarity_transformation.scale,
              similarity_transformation.scale,
              10 * tolerance);
}

TEST(EstimateSimilarityTransformation2D3D, AllInliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 1.0;
  const double kNoise = 0.0;
  const double kPoseTolerance = 1e-4;

  SimilarityTransformation similarity_transformation;
  similarity_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                                            .normalized()).toRotationMatrix();
  similarity_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  similarity_transformation.scale = 0.8;
  ExecuteRandomTest(options,
                    similarity_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance);
}

TEST(EstimateSimilarityTransformation2D3D, AllInliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 1.0;
  const double kNoise = 1.0;
  const double kPoseTolerance = 1e-2;

  SimilarityTransformation similarity_transformation;
  similarity_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                                            .normalized()).toRotationMatrix();
  similarity_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  similarity_transformation.scale = 0.8;
  ExecuteRandomTest(options,
                    similarity_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance);
}

TEST(EstimateSimilarityTransformation2D3D, OutliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 0.8;
  const double kNoise = 0.0;
  const double kPoseTolerance = 1e-4;

  SimilarityTransformation similarity_transformation;
  similarity_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                                            .normalized()).toRotationMatrix();
  similarity_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  similarity_transformation.scale = 0.8;
  ExecuteRandomTest(options,
                    similarity_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance);
}

TEST(EstimateSimilarityTransformation2D3D, OutliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 0.8;
  const double kNoise = 1.0;
  const double kPoseTolerance = 1e-2;

  SimilarityTransformation similarity_transformation;
  similarity_transformation.rotation =
      Eigen::AngleAxisd(DegToRad(12.0), Eigen::Vector3d(1.0, 0.2, -0.8)
                                            .normalized()).toRotationMatrix();
  similarity_transformation.translation = Eigen::Vector3d(-1.3, 2.1, 0.5);
  similarity_transformation.scale = 0.8;
  ExecuteRandomTest(options,
                    similarity_transformation,
                    kInlierRatio,
                    kNoise,
                    kPoseTolerance);
}

}  // namespace theia
