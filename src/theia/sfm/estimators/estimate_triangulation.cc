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

#include "theia/sfm/estimators/estimate_triangulation.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/sfm/types.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/ransac.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {

namespace {
// The pixel observation and projection matrix needed in order to triangulate a
// 3D point.
struct PointObservation {
  Matrix3x4d projection_matrix;
  Camera camera;
  Eigen::Vector2d normalized_feature;
  Eigen::Vector2d observed_pixel;
};

// Returns true if the point is in front of the camera and false if the point is
// behind the camera.
bool IsPointInFrontOfCamera(const Matrix3x4d& projection_matrix,
                            const Eigen::Vector4d& point) {
  return point.dot(projection_matrix.row(2)) > 0;
}

class TriangulationEstimator
    : public Estimator<PointObservation, Eigen::Vector4d> {
 public:
  TriangulationEstimator() {}

  double SampleSize() const { return 2; }

  // Triangulates the 3D point from 2 observations.
  bool EstimateModel(const std::vector<PointObservation>& observations,
                     std::vector<Eigen::Vector4d>* triangulated_points) const {
    triangulated_points->resize(1);
    if (!Triangulate(observations[0].projection_matrix,
                     observations[1].projection_matrix,
                     observations[0].normalized_feature,
                     observations[1].normalized_feature,
                     &triangulated_points->at(0))) {
      return false;
    }

    // Only return true if the point is in front of both cameras and the
    // triangulation was a success.
    return IsPointInFrontOfCamera(observations[0].projection_matrix,
                                  triangulated_points->at(0)) &&
           IsPointInFrontOfCamera(observations[1].projection_matrix,
                                  triangulated_points->at(0));
  }

  double Error(const PointObservation& observation,
               const Eigen::Vector4d& triangulated_point) const {
    Eigen::Vector2d reprojection;
    const double depth =
        observation.camera.ProjectPoint(triangulated_point, &reprojection);
    if (depth <= 0) {
      return std::numeric_limits<double>::max();
    }
    return (observation.observed_pixel - reprojection).squaredNorm();
  }
};

}  // namespace

bool EstimateTriangulation(const RansacParameters& ransac_params,
                           const std::vector<Camera>& cameras,
                           const std::vector<Eigen::Vector2d>& features,
                           Eigen::Vector4d* triangulated_point,
                           RansacSummary* summary) {
  CHECK_EQ(cameras.size(), features.size());
  CHECK_NOTNULL(triangulated_point);

  // If we only have a few data points, then we should exhaustively search the
  // solution space for the best combination.
  static const int kMaxNumDataPointsForExhaustiveSearch = 15;
  if (cameras.size() < 2) {
    return false;
  }

  // Create point correspondences.
  std::vector<PointObservation> point_observations(cameras.size());
  for (int i = 0; i < point_observations.size(); i++) {
    // Create the projection atirx.
    Matrix3x4d projection_matrix;
    projection_matrix.leftCols<3>() =
        cameras[i].GetOrientationAsRotationMatrix();
    projection_matrix.rightCols<1>() =
        -projection_matrix.leftCols<3>() * cameras[i].GetPosition();

    point_observations[i].projection_matrix = projection_matrix;
    point_observations[i].camera = cameras[i];
    point_observations[i].normalized_feature =
        cameras[i].PixelToNormalizedCoordinates(features[i]).hnormalized();
    point_observations[i].observed_pixel = features[i];
  }

  // RANSAC triangulation.
  TriangulationEstimator triangulation_estimator;
  std::unique_ptr<SampleConsensusEstimator<TriangulationEstimator> > ransac;
  if (cameras.size() <= kMaxNumDataPointsForExhaustiveSearch) {
    // Set the minimum number of iterations to be the exact number of possible
    // combinations. This forces all combinations to be tested.
    const int num_combinations = cameras.size() * (cameras.size() - 1) / 2;
    RansacParameters exhaustive_params = ransac_params;
    exhaustive_params.min_iterations = num_combinations;
    exhaustive_params.max_iterations = num_combinations;
    ransac = CreateAndInitializeRansacVariant<TriangulationEstimator>(
        RansacType::EXHAUSTIVE, exhaustive_params, triangulation_estimator);
  } else {
    ransac = CreateAndInitializeRansacVariant<TriangulationEstimator>(
        RansacType::RANSAC, ransac_params, triangulation_estimator);
  }

  // Run the RANSAC scheme.
  CHECK(ransac->Initialize());
  return ransac->Estimate(point_observations, triangulated_point, summary);
}

}  // namespace theia
