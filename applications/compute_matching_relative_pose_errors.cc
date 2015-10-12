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
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>

#include <memory>
#include <string>

DEFINE_string(
    matches, "",
    "Matches file that has been written with WriteMatchesAndGeometry");

DEFINE_string(reconstruction,
              "",
              "Reconstruction to use as ground truth.");

using theia::Reconstruction;
using theia::TrackId;
using theia::ViewId;

// Computes the error in the relative rotation based on the ground truth
// rotations rotation1 and rotation2 (which specify world-to-camera
// transformations).
double ComputeRelativeRotationError(const Eigen::Matrix3d& rotation1,
                                    const Eigen::Matrix3d& rotation2,
                                    const Eigen::Vector3d& relative_rotation) {
  Eigen::Matrix3d relative_rotation_matrix;
  ceres::AngleAxisToRotationMatrix(
      relative_rotation.data(),
      ceres::ColumnMajorAdapter3x3(relative_rotation_matrix.data()));
  const Eigen::Matrix3d loop_rotation = relative_rotation_matrix.transpose() *
                                        (rotation2 * rotation1.transpose());
  Eigen::Vector3d loop_rotation_aa;
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(loop_rotation.data()),
      loop_rotation_aa.data());
  return theia::RadToDeg(loop_rotation_aa.norm());
}

double ComputeRelativeTranslationError(
    const Eigen::Vector3d& position1,
    const Eigen::Vector3d& position2,
    const Eigen::Matrix3d& rotation1,
    const Eigen::Vector3d& relative_translation) {
  const Eigen::Vector3d world_translation =
      rotation1 * (position2 - position1).normalized();
  return theia::RadToDeg(acos(
      theia::Clamp(relative_translation.dot(world_translation), -1.0, 1.0)));
}

void EvaluateRelativeError(
    const std::vector<std::string>& view_names,
    const std::vector<theia::ImagePairMatch>& matches,
    const Reconstruction& gt_reconstruction) {
  // For each edge, get the rotate translation and check the error.
  int num_matches_evaluated = 0;
  const std::vector<double> histogram_bins = {2,  5,   10,  15,  25,  50,
                                              90, 135, 180, 225, 270, 316};

  theia::PoseError relative_pose_error(histogram_bins, histogram_bins);
  for (const auto& match : matches) {
    const std::string view1_name = match.image1;
    const std::string view2_name = match.image2;

    const ViewId view_id1 =
        gt_reconstruction.ViewIdFromName(view1_name);
    const ViewId view_id2 =
        gt_reconstruction.ViewIdFromName(view2_name);
    const theia::View* view1 = gt_reconstruction.View(view_id1);
    const theia::View* view2 = gt_reconstruction.View(view_id2);
    if (view1 == nullptr || view2 == nullptr) {
      continue;
    }
    const theia::Camera& camera1 = view1->Camera();
    const theia::Camera& camera2 = view2->Camera();

    const double rotation_angular_error =
        ComputeRelativeRotationError(camera1.GetOrientationAsRotationMatrix(),
                                     camera2.GetOrientationAsRotationMatrix(),
                                     match.twoview_info.rotation_2);
    const double translation_angular_error = ComputeRelativeTranslationError(
        camera1.GetPosition(),
        camera2.GetPosition(),
        camera1.GetOrientationAsRotationMatrix(),
        match.twoview_info.position_2);

    relative_pose_error.AddError(rotation_angular_error,
                                 translation_angular_error);
    ++num_matches_evaluated;
  }
  LOG(INFO) << "Evaluated " << view_names.size() << " common views containing "
            << num_matches_evaluated << " two-view matches.";
  LOG(INFO) << relative_pose_error.PrintMeanMedianHistogram();
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  LOG(INFO) << "Reading the matches.";

  std::vector<std::string> view_names;
  std::vector<theia::CameraIntrinsicsPrior> camera_intrinsics_prior;
  std::vector<theia::ImagePairMatch> matches;
  CHECK(theia::ReadMatchesAndGeometry(FLAGS_matches,
                                      &view_names,
                                      &camera_intrinsics_prior,
                                      &matches))
      << "Could not read matches from " << FLAGS_matches;

  std::unique_ptr<theia::Reconstruction> reconstruction(
      new theia::Reconstruction());
  CHECK(theia::ReadReconstruction(FLAGS_reconstruction, reconstruction.get()))
      << "Could not read reconstruction from " << FLAGS_reconstruction;

  EvaluateRelativeError(view_names, matches, *reconstruction);

  return 0;
}
