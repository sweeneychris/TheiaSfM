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

#include "theia/sfm/transformation/align_reconstructions.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <string>
#include <vector>

#include "theia/sfm/find_common_views_by_name.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/transformation/align_point_clouds.h"
#include "theia/sfm/transformation/transform_reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/solvers/ransac.h"
#include "theia/solvers/sample_consensus_estimator.h"

namespace theia {
namespace {

struct CameraCorrespondence {
  Eigen::Vector3d camera1;
  Eigen::Vector3d camera2;
};

struct SimilarityTransformation {
  Eigen::Matrix3d rotation;
  Eigen::Vector3d translation;
  double scale;
};

class CameraAlignmentEstimator
    : public Estimator<CameraCorrespondence, SimilarityTransformation> {
 public:
  CameraAlignmentEstimator() {}

  double SampleSize() const { return 4; }

  bool EstimateModel(
      const std::vector<CameraCorrespondence>& correspondences,
      std::vector<SimilarityTransformation>* sim_transforms) const {
    SimilarityTransformation sim_transform;
    std::vector<Eigen::Vector3d> positions1(correspondences.size());
    std::vector<Eigen::Vector3d> positions2(correspondences.size());
    for (int i = 0; i < correspondences.size(); i++) {
      positions1[i] = correspondences[i].camera1;
      positions2[i] = correspondences[i].camera2;
    }

    AlignPointCloudsUmeyama(positions2,
                            positions1,
                            &sim_transform.rotation,
                            &sim_transform.translation,
                            &sim_transform.scale);
    sim_transforms->emplace_back(sim_transform);
    return true;
  }

  double Error(const CameraCorrespondence& cameras,
               const SimilarityTransformation& sim_transform) const {
    const Eigen::Vector3d transformed_camera =
        sim_transform.scale * sim_transform.rotation * cameras.camera2 +
        sim_transform.translation;
    return (cameras.camera1 - transformed_camera).squaredNorm();
  }
};

}  // namespace

void AlignReconstructions(const Reconstruction& reconstruction1,
                          Reconstruction* reconstruction2) {
  CHECK_NOTNULL(reconstruction2);

  const std::vector<std::string> common_view_names =
      FindCommonViewsByName(reconstruction1, *reconstruction2);

  // Collect the positions of all common views.
  std::vector<Eigen::Vector3d> positions1(common_view_names.size());
  std::vector<Eigen::Vector3d> positions2(common_view_names.size());
  for (int i = 0; i < common_view_names.size(); i++) {
    const ViewId view_id1 =
        reconstruction1.ViewIdFromName(common_view_names[i]);
    const ViewId view_id2 =
        reconstruction2->ViewIdFromName(common_view_names[i]);
    positions1[i] =
        reconstruction1.View(view_id1)->Camera().GetPosition();
    positions2[i] =
        reconstruction2->View(view_id2)->Camera().GetPosition();
  }

  // Align the positions.
  Eigen::Matrix3d rotation;
  Eigen::Vector3d translation;
  double scale;
  AlignPointCloudsUmeyama(positions2,
                          positions1,
                          &rotation,
                          &translation,
                          &scale);

  // Apply the similarity transformation to the reconstruction.
  TransformReconstruction(rotation, translation, scale, reconstruction2);
}

void AlignReconstructionsRobust(
    const double robust_error_threshold,
    const Reconstruction& reconstruction1,
    Reconstruction* reconstruction2) {
  CHECK_NOTNULL(reconstruction2);

  const std::vector<std::string> common_view_names =
      FindCommonViewsByName(reconstruction1, *reconstruction2);

  // Collect the positions of all common views.
  std::vector<CameraCorrespondence> correspondences(common_view_names.size());
  for (int i = 0; i < common_view_names.size(); i++) {
    const ViewId view_id1 =
        reconstruction1.ViewIdFromName(common_view_names[i]);
    const ViewId view_id2 =
        reconstruction2->ViewIdFromName(common_view_names[i]);
    correspondences[i].camera1 =
        reconstruction1.View(view_id1)->Camera().GetPosition();
    correspondences[i].camera2 =
        reconstruction2->View(view_id2)->Camera().GetPosition();
  }

  // Estimate with RANSAC.
  RansacParameters params;
  params.max_iterations = 1000;
  params.use_mle = true;
  params.error_thresh = robust_error_threshold * robust_error_threshold;
  params.failure_probability = 1e-4;

  CameraAlignmentEstimator estimator;
  Ransac<CameraAlignmentEstimator> ransac(params, estimator);
  CHECK(ransac.Initialize()) << "Could not initialize RANSAC for similarity "
                                "transformation estimation.";
  SimilarityTransformation sim_transform;
  RansacSummary summary;
  CHECK(ransac.Estimate(correspondences, &sim_transform, &summary))
      << "Could not align models with RANSAC";
  CHECK_GT(summary.inliers.size(), 0)
      << "No inliers could be found for estimating the similarity "
         "transformation. Try using a higher error threshold.";

  // Align the reconstructions using the inliers.
  std::vector<Eigen::Vector3d> positions1(summary.inliers.size());
  std::vector<Eigen::Vector3d> positions2(summary.inliers.size());
  for (int i = 0; i < summary.inliers.size(); i++) {
    positions1[i] = correspondences[summary.inliers[i]].camera1;
    positions2[i] = correspondences[summary.inliers[i]].camera2;
  }

  // Align the positions.
  Eigen::Matrix3d rotation;
  Eigen::Vector3d translation;
  double scale;
  AlignPointCloudsUmeyama(positions2,
                          positions1,
                          &rotation,
                          &translation,
                          &scale);

  // Apply the similarity transformation to the reconstruction.
  TransformReconstruction(rotation, translation, scale, reconstruction2);
}

}  // namespace theia
