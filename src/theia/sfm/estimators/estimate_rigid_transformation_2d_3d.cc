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

#include "theia/sfm/estimators/estimate_rigid_transformation_2d_3d.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
#include <vector>
#include <glog/logging.h>

#include "theia/sfm/estimators/camera_and_feature_correspondence_2d_3d.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/pose/upnp.h"
#include "theia/sfm/rigid_transformation.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {
namespace {

// An estimator for computing the absolute pose from 2 feature
// correspondences. The feature correspondences should be normalized by the
// focal length with the principal point at (0, 0) and rotated into the same
// coordinate system as the points.
class NonCentralCameraPoseEstimator
    : public Estimator<CameraAndFeatureCorrespondence2D3D,
                       RigidTransformation> {
 public:
  NonCentralCameraPoseEstimator() {
    estimator_.reset(new class Upnp);
  }
  ~NonCentralCameraPoseEstimator() = default;

  // Get the minimum number of samples needed to generate a model.
  virtual double SampleSize() const override { return 4; }

  // Given a set of data points, estimate the model. Users should implement this
  // function appropriately for the task being solved. Returns true for
  // successful model estimation (and outputs model), false for failed
  // estimation. Typically, this is a minimal set, but it is not required to be.
  bool EstimateModel(
      const std::vector<CameraAndFeatureCorrespondence2D3D>& data,
      std::vector<RigidTransformation>* models) const override {
    CHECK_GE(data.size(), static_cast<size_t>(SampleSize()));
    CHECK_NOTNULL(models)->clear();
    // Input datum for Upnp.
    std::vector<Eigen::Vector3d> ray_directions;
    std::vector<Eigen::Vector3d> ray_origins;
    std::vector<Eigen::Vector3d> world_points;

    // Prepare data.
    ray_directions.reserve(data.size());
    ray_origins.reserve(data.size());
    world_points.reserve(data.size());
    for (const CameraAndFeatureCorrespondence2D3D& correspondence : data) {
      ray_directions.emplace_back(correspondence.camera.PixelToUnitDepthRay(
          correspondence.observation).normalized());
      ray_origins.emplace_back(correspondence.camera.GetPosition());
      world_points.emplace_back(correspondence.point3d.hnormalized());
    }

    // Estimate pose.
    std::vector<Eigen::Quaterniond> soln_rotations;
    std::vector<Eigen::Vector3d> soln_translations;
    if (!estimator_->EstimatePose(ray_origins,
                                  ray_directions,
                                  world_points,
                                  &soln_rotations,
                                  &soln_translations)) {
      return false;
    }

    // Generate models.
    models->resize(soln_rotations.size());
    for (int i = 0; i < soln_rotations.size(); ++i) {
      models->at(i).rotation = soln_rotations[i].toRotationMatrix();
      models->at(i).translation = soln_translations[i];
    }
    return true;
  }

  // Given a model and a data point, calculate the error. Users should implement
  // this function appropriately for the task being solved.
  double Error(const CameraAndFeatureCorrespondence2D3D& data,
               const RigidTransformation& model) const override {
    Eigen::Vector2d reprojection;
    // If the point is reprojected behind the camera, return the maximum
    // possible error.
    const Eigen::Vector4d transformed_point =
        (model.rotation * data.point3d.hnormalized() +
         model.translation).homogeneous();
    if (data.camera.ProjectPoint(transformed_point, &reprojection) < 0) {
      return std::numeric_limits<double>::max();
    }

    // Return the squared reprojection error.
    return (data.observation - reprojection).squaredNorm();
  }
  
 private:
  // Upnp estimator.
  std::unique_ptr<class Upnp> estimator_;
  
  DISALLOW_COPY_AND_ASSIGN(NonCentralCameraPoseEstimator);
};

}  // namespace

bool EstimateRigidTransformation2D3D(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<CameraAndFeatureCorrespondence2D3D>& correspondences,
    RigidTransformation* estimated_pose,
    RansacSummary* ransac_summary) {
  NonCentralCameraPoseEstimator pose_estimator;
  std::unique_ptr<SampleConsensusEstimator<NonCentralCameraPoseEstimator>>
      ransac = CreateAndInitializeRansacVariant(ransac_type,
                                                ransac_params,
                                                pose_estimator);
  // Estimate the pose.
  return ransac->Estimate(correspondences,
                          estimated_pose,
                          ransac_summary);
}

bool EstimateRigidTransformation2D3D(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence2D3D>& normalized_correspondences,
    RigidTransformation* estimated_pose,
    RansacSummary* ransac_summary) {
  // const int num_correspondences = normalized_correspondences.size();
  // std::vector<NonCentralCameraFeatureCorrespondence> correspondences;
  // correspondences.reserve(num_correspondences);
  // for (int i = 0; i < num_correspondences; ++i) {
  //   const FeatureCorrespondence2D3D& normalized_correspondence =
  //       normalized_correspondences[i];
  //   correspondences.emplace_back(
  //       normalized_correspondence.feature.homogeneous().normalized(),
  //       Eigen::Vector3d::Zero(),
  //       normalized_correspondence.world_point);
  // }
  // return EstimateRigidTransformation2D3D(ransac_params,
  //                                        ransac_type,
  //                                        correspondences,
  //                                        estimated_pose,
  //                                        ransac_summary);
  return false;
}

}  // namespace theia
