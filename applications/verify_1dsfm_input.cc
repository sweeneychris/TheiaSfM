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

DEFINE_string(1dsfm_dataset_directory, "",
              "Dataset where the 1dSFM dataset is located. Do not include a "
              "trailing slash.");

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
    const theia::ViewGraph& view_graph,
    const Reconstruction& reconstruction_1dsfm,
    const Reconstruction& gt_reconstruction) {
  // For each edge, get the rotate translation and check the error.
  std::vector<double> histogram_bins = {2,  5,   10,  15,  25,  50,
                                        90, 135, 180, 225, 270, 315};
  theia::PoseError pose_error(histogram_bins, histogram_bins);

  const auto& edges = view_graph.GetAllEdges();
  for (const auto& edge : edges) {
    // The reconstruction/view graph from 1dSFM may have a different mapping of
    // names to ViewIds than the ground truth reconstruction, so we have to do a
    // name lookup in the ground truth reconstruction to ensure we have the same
    // view.
    const theia::View* view1_1dsfm =
        reconstruction_1dsfm.View(edge.first.first);
    const theia::View* view2_1dsfm =
        reconstruction_1dsfm.View(edge.first.second);
    if (view1_1dsfm == nullptr || view2_1dsfm == nullptr) {
      continue;
    }

    const ViewId view_id1 =
        gt_reconstruction.ViewIdFromName(view1_1dsfm->Name());
    const ViewId view_id2 =
        gt_reconstruction.ViewIdFromName(view2_1dsfm->Name());
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
                                     edge.second.rotation_2);

    const double translation_angular_error = ComputeRelativeTranslationError(
        camera1.GetPosition(),
        camera2.GetPosition(),
        camera1.GetOrientationAsRotationMatrix(),
        edge.second.position_2);
    pose_error.AddError(rotation_angular_error, translation_angular_error);
  }

  LOG(INFO) << "Relative pose errors for 1dsfm = \n"
            << pose_error.PrintMeanMedianHistogram();
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  LOG(INFO) << "Reading the relative poses from the 1dsfm dataset.";
  std::unique_ptr<Reconstruction> reconstruction_1dsfm(new Reconstruction);
  std::unique_ptr<theia::ViewGraph> view_graph_1dsfm(new theia::ViewGraph);
  CHECK(Read1DSFM(FLAGS_1dsfm_dataset_directory,
                  reconstruction_1dsfm.get(),
                  view_graph_1dsfm.get()))
      << "Could not read 1dsfm dataset from " << FLAGS_1dsfm_dataset_directory;

  const std::string lists_file = FLAGS_1dsfm_dataset_directory + "/list.txt";
  const std::string bundle_file =
      FLAGS_1dsfm_dataset_directory + "/gt_bundle.out";
  std::unique_ptr<theia::Reconstruction> gt_reconstruction(
      new theia::Reconstruction());
  LOG(INFO) << "Converting ground truth bundler file to Theia reconstruction.";
  CHECK(
      theia::ReadBundlerFiles(lists_file, bundle_file, gt_reconstruction.get()))
      << "Could not the ground truth Bundler file at " << bundle_file;

  EvaluateRelativeError(*view_graph_1dsfm,
                        *reconstruction_1dsfm,
                        *gt_reconstruction);

  return 0;
}
