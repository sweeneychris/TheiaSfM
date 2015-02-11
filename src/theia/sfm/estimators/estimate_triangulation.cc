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

#include "theia/solvers/estimator.h"
#include "theia/solvers/ransac.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/sfm/types.h"

namespace theia {

namespace {
// The pixel observation and projection matrix needed in order to triangulate a
// 3D point.
struct PointObservation {
  Matrix3x4d projection_matrix;
  Eigen::Vector2d feature;
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
    // TODO(cmsweeney): We do not check the angle between the two views at the
    // moment. This requires the ray direction of each feature meaning we would
    // have to either decompose the projection matrix or pass in the ray
    // direction as part of the Point Observation. RANSAC should be good enough
    // at filtering out these bad solutions so we ignore this for now.
    triangulated_points->resize(1);
    if (!Triangulate(observations[0].projection_matrix,
                     observations[1].projection_matrix,
                     observations[0].feature,
                     observations[1].feature,
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
    if (!IsPointInFrontOfCamera(observation.projection_matrix,
                                triangulated_point)) {
      return std::numeric_limits<double>::max();
    }

    const Eigen::Vector2d reprojection =
        (observation.projection_matrix * triangulated_point).hnormalized();
    return (observation.feature - reprojection).squaredNorm();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(TriangulationEstimator);
};

}  // namespace

bool EstimateTriangulation(const RansacParameters& ransac_params,
                           const std::vector<Matrix3x4d>& projection_matrices,
                           const std::vector<Eigen::Vector2d>& features,
                           Eigen::Vector4d* triangulated_point,
                           RansacSummary* summary) {
  if (projection_matrices.size() < 2) {
    return false;
  }

  // Create point correspondences.
  std::vector<PointObservation> point_observations(
      projection_matrices.size());
  for (int i = 0; i < point_observations.size(); i++) {
    point_observations[i].projection_matrix = projection_matrices[i];
    point_observations[i].feature = features[i];
  }

  // RANSAC triangulation.
  TriangulationEstimator triangulation_estimator;
  Ransac<TriangulationEstimator> ransac(ransac_params,
                                        triangulation_estimator);
  CHECK(ransac.Initialize());
  if (!ransac.Estimate(point_observations, triangulated_point, summary)) {
    return false;
  }

  return true;
}

}  // namespace theia
