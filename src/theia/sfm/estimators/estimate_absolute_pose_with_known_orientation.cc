#include "theia/sfm/estimators/estimate_absolute_pose_with_known_orientation.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/pose/position_from_two_rays.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {
namespace {

// Rotate the correspondences to be in the world space instead of camera
// coordinates.
void RotateCorrespondences(
    const std::vector<FeatureCorrespondence2D3D>& normalized_correspondences,
    const Eigen::Vector3d& camera_orientation,
    std::vector<FeatureCorrespondence2D3D>* rotated_correspondences) {
  Eigen::Matrix3d camera_to_world_rotation;
  ceres::AngleAxisToRotationMatrix(
      camera_orientation.data(),
      ceres::ColumnMajorAdapter3x3(camera_to_world_rotation.data()));
  camera_to_world_rotation.transposeInPlace();

  rotated_correspondences->resize(normalized_correspondences.size());
  for (int i = 0; i < normalized_correspondences.size(); i++) {
    const Eigen::Vector3d feature_ray =
        camera_to_world_rotation *
        normalized_correspondences[i].feature.homogeneous();
    (*rotated_correspondences)[i].feature = feature_ray.hnormalized();
    (*rotated_correspondences)[i].world_point =
        normalized_correspondences[i].world_point;
  }
}

// An estimator for computing the absolute pose from 2 feature
// correspondences. The feature correspondences should be normalized by the
// focal length with the principal point at (0, 0) and rotated into the same
// coordinate system as the points.
class AbsolutePoseWithKnownOrientationEstimator
    : public Estimator<FeatureCorrespondence2D3D, Eigen::Vector3d> {
 public:
  AbsolutePoseWithKnownOrientationEstimator() {}

  // 2 correspondences are needed to determine the absolute position.
  double SampleSize() const { return 2; }

  // Estimates candidate absolute poses from correspondences.
  bool EstimateModel(
      const std::vector<FeatureCorrespondence2D3D>& correspondences,
      std::vector<Eigen::Vector3d>* absolute_positions) const {
    Eigen::Vector3d position;
    if (!PositionFromTwoRays(correspondences[0].feature,
                             correspondences[0].world_point,
                             correspondences[1].feature,
                             correspondences[1].world_point,
                             &position)) {
      return false;
    }

    absolute_positions->emplace_back(position);
    return true;
  }

  // The error for a correspondences given an absolute position. This is the
  // squared reprojection error.
  double Error(const FeatureCorrespondence2D3D& correspondence,
               const Eigen::Vector3d& absolute_position) const {
    // The reprojected point is computed as X - c where c is the position, and X
    // is the 3D point.
    const Eigen::Vector2d reprojected_feature =
        (correspondence.world_point - absolute_position).hnormalized();
    return (reprojected_feature - correspondence.feature).squaredNorm();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(AbsolutePoseWithKnownOrientationEstimator);
};

}  // namespace

// Estimates the absolute pose with known orientation. It is assumed that the 2D
// features in the 2D-3D correspondences are normalized by the camera
// intrinsics. Returns true if the position could be successfully estimated and
// false otherwise. The quality of the result depends on the quality of the
// input data.
bool EstimateAbsolutePoseWithKnownOrientation(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const Eigen::Vector3d& camera_orientation,
    const std::vector<FeatureCorrespondence2D3D>& normalized_correspondences,
    Eigen::Vector3d* camera_position,
    RansacSummary* ransac_summary) {

  // Rotate the correspondences.
  std::vector<FeatureCorrespondence2D3D> rotated_correspondences;
  RotateCorrespondences(normalized_correspondences,
                        camera_orientation,
                        &rotated_correspondences);

  AbsolutePoseWithKnownOrientationEstimator absolute_pose_estimator;
  std::unique_ptr<
      SampleConsensusEstimator<AbsolutePoseWithKnownOrientationEstimator> >
      ransac = CreateAndInitializeRansacVariant(ransac_type,
                                                ransac_params,
                                                absolute_pose_estimator);
  // Estimate the absolute pose.
  return ransac->Estimate(rotated_correspondences,
                          camera_position,
                          ransac_summary);
}

}  // namespace theia
