#include "theia/sfm/estimators/estimate_absolute_pose_with_known_orientation.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/pose/essential_matrix_utils.h"
#include "theia/sfm/pose/util.h"
#include "theia/sfm/pose/relative_pose_from_two_points_with_known_rotation.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/util/util.h"

namespace theia {
namespace {

// An estimator for computing the relative pose from 2 feature
// correspondences. The feature correspondences should be normalized by the
// focal length with the principal point at (0, 0).
class RelativePoseWithKnownOrientationEstimator
    : public Estimator<FeatureCorrespondence, Eigen::Vector3d> {
 public:
  RelativePoseWithKnownOrientationEstimator() {}

  // 2 correspondences are needed to determine the relative position.
  double SampleSize() const { return 2; }

  // Estimates candidate relative poses from correspondences.
  bool EstimateModel(const std::vector<FeatureCorrespondence>& correspondences,
                     std::vector<Eigen::Vector3d>* relative_positions) const {
    Eigen::Vector3d position;
    const Eigen::Vector2d rotated_features1[2] = {correspondences[0].feature1,
                                                  correspondences[1].feature1};
    const Eigen::Vector2d rotated_features2[2] = {correspondences[0].feature2,
                                                  correspondences[1].feature2};
    if (!RelativePoseFromTwoPointsWithKnownRotation(
            rotated_features1, rotated_features2, &position)) {
      return false;
    }

    relative_positions->emplace_back(position);
    return true;
  }

  // The error for a correspondences given an relative position. This is the
  // squared reprojection error.
  double Error(const FeatureCorrespondence& correspondence,
               const Eigen::Vector3d& relative_position) const {
    static const Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    return SquaredSampsonDistance(CrossProductMatrix(-relative_position),
                                  correspondence.feature1,
                                  correspondence.feature2);
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(RelativePoseWithKnownOrientationEstimator);
};

}  // namespace

bool EstimateRelativePoseWithKnownOrientation(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence>& rotated_correspondences,
    Eigen::Vector3d* relative_camera2_position,
    RansacSummary* ransac_summary) {
  RelativePoseWithKnownOrientationEstimator relative_pose_estimator;
  std::unique_ptr<
      SampleConsensusEstimator<RelativePoseWithKnownOrientationEstimator> >
      ransac = CreateAndInitializeRansacVariant(ransac_type,
                                                ransac_params,
                                                relative_pose_estimator);
  // Estimate the relative pose.
  return ransac->Estimate(rotated_correspondences,
                          relative_camera2_position,
                          ransac_summary);
}

}  // namespace theia
