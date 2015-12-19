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

#include "theia/sfm/estimators/estimate_similarity_transformation_2d_3d.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <limits>
#include <memory>
#include <vector>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/similarity_transformation.h"
#include "theia/sfm/transformation/gdls_similarity_transform.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {
namespace {

inline void TransformCamera(const SimilarityTransformation& sim_transform,
                            Camera* camera) {
  const Eigen::Vector3d old_position = camera->GetPosition();
  const Eigen::Vector3d new_position =
      sim_transform.scale * sim_transform.rotation *
      old_position + sim_transform.translation;
  camera->SetPosition(new_position);

  const Eigen::Matrix3d old_orientation =
      camera->GetOrientationAsRotationMatrix();
  const Eigen::Matrix3d new_orientation =
      old_orientation * sim_transform.rotation.transpose();
  camera->SetOrientationFromRotationMatrix(new_orientation);
}

// An estimator for computing the similarity transformation from 4 2D-3D
// correspondences using the gDLS algorithm to minimize reprojection error.
class GdlsSimilarityTransformationEstimator
    : public Estimator<CameraAndFeatureCorrespondence2D3D,
                       SimilarityTransformation> {
 public:
  GdlsSimilarityTransformationEstimator() {}

  // 3 correspondences are needed to determine the absolute pose.
  double SampleSize() const { return 4; }

  // Estimates candidate absolute poses from correspondences.
  bool EstimateModel(
      const std::vector<CameraAndFeatureCorrespondence2D3D>& correspondences,
      std::vector<SimilarityTransformation>* similarity_transformations) const {
    std::vector<Eigen::Vector3d> ray_origins(4), ray_directions(4),
        world_points(4);
    for (int i = 0; i < 4; i++) {
      ray_origins[i] = correspondences[i].camera.GetPosition();
      ray_directions[i] = correspondences[i].camera.PixelToUnitDepthRay(
          correspondences[i].observation).normalized();
      world_points[i] = correspondences[i].point3d.hnormalized();
    }

    // Compute the similarity transformation. Note that this function computes
    // R, t, and s such that:
    //
    //   s * c_i + alpha_i * x_i = R * X_i + t
    //
    // where c_i is the camera position, alpha_i is the depth of the feature,
    // x_i is the unit-norm feature observation, and X_i is the 3D point.
    std::vector<Eigen::Quaterniond> rotations;
    std::vector<Eigen::Vector3d> translations;
    std::vector<double> scales;
    GdlsSimilarityTransform(ray_origins,
                            ray_directions,
                            world_points,
                            &rotations,
                            &translations,
                            &scales);

    // Aggregate the solutions, modifying the output so that R, t, s are of the
    // more useful form of:
    //
    //   s * R * (c_i + alpha_i * x_i) + t = X_i
    //
    // which transforms only the camera coordinate system so that it is aligned
    // with the 3D points.
    for (int i = 0; i < rotations.size(); i++) {
      SimilarityTransformation similarity_transformation;
      similarity_transformation.rotation =
          rotations[i].toRotationMatrix().transpose();
      similarity_transformation.translation =
          similarity_transformation.rotation * -translations[i];
      similarity_transformation.scale = scales[i];
      similarity_transformations->emplace_back(similarity_transformation);
    }
    return similarity_transformations->size() > 0;
  }

  // The error for a correspondences given an absolute pose. This is the squared
  // reprojection error.
  double Error(
      const CameraAndFeatureCorrespondence2D3D& correspondence,
      const SimilarityTransformation& similarity_transformation) const {
    // Apply the similarity transformation to the camera.
    Camera transformed_camera = correspondence.camera;
    TransformCamera(similarity_transformation, &transformed_camera);

    Eigen::Vector2d reprojection;
    // If the point is reprojected behind the camera, return the maximum
    // possible error.
    if (transformed_camera.ProjectPoint(correspondence.point3d, &reprojection) <
        0) {
      return std::numeric_limits<double>::max();
    }

    // Return the squared reprojection error.
    return (correspondence.observation - reprojection).squaredNorm();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(GdlsSimilarityTransformationEstimator);
};

}  // namespace

bool EstimateSimilarityTransformation2D3D(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<CameraAndFeatureCorrespondence2D3D>& correspondences,
    SimilarityTransformation* similarity_transformation,
    RansacSummary* ransac_summary) {
  GdlsSimilarityTransformationEstimator similarity_transformation_estimator;
  std::unique_ptr <
      SampleConsensusEstimator<GdlsSimilarityTransformationEstimator> > ransac =
      CreateAndInitializeRansacVariant(ransac_type,
                                       ransac_params,
                                       similarity_transformation_estimator);
  // Estimate the absolute pose.
  return ransac->Estimate(correspondences,
                          similarity_transformation,
                          ransac_summary);
}

}  // namespace theia
