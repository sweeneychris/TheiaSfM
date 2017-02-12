// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "theia/io/read_matches.h"
#include "theia/io/reconstruction_reader.h"
#include "theia/matching/image_pair_match.h"
#include "theia/sfm/find_common_views_by_name.h"
#include "theia/sfm/incremental_reconstruction_estimator.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/transformation/align_reconstructions.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {
void ReadInput(Reconstruction* gt_reconstruction,
               Reconstruction* reconstruction,
               ViewGraph* view_graph) {
  const std::string gt_reconstruction_filename =
      THEIA_DATA_DIR + std::string("/sfm/gt_fountain11.bin");
  const std::string reconstruction_filename =
      THEIA_DATA_DIR + std::string("/sfm/fountain11.bin");
  const std::string matches_filename =
      THEIA_DATA_DIR + std::string("/sfm/fountain11_matches.bin");

  // Read the reconstruction file.
  CHECK(ReadReconstruction(gt_reconstruction_filename, gt_reconstruction));
  CHECK(ReadReconstruction(reconstruction_filename, reconstruction));

  // Set the views and tracks to be unestimated so that we can re-estimate them
  // with incremental SfM.
  const auto view_ids = reconstruction->ViewIds();
  for (const ViewId view_id : view_ids) {
    reconstruction->MutableView(view_id)->SetEstimated(false);
  }
  const auto track_ids = reconstruction->TrackIds();
  for (const TrackId track_id : track_ids) {
    reconstruction->MutableTrack(track_id)->SetEstimated(false);
  }

  // Read the matches file.
  std::vector<std::string> image_files;
  std::vector<CameraIntrinsicsPrior> camera_intrinsics_prior;
  std::vector<ImagePairMatch> image_matches;

  // Read in match file.
  CHECK(ReadMatchesAndGeometry(matches_filename,
                               &image_files,
                               &camera_intrinsics_prior,
                               &image_matches));

  // Add the matches to the view graph.
  for (const ImagePairMatch& match : image_matches) {
    TwoViewInfo info = match.twoview_info;
    const ViewId view_id1 = reconstruction->ViewIdFromName(match.image1);
    const ViewId view_id2 = reconstruction->ViewIdFromName(match.image2);
    if (view_id1 == kInvalidViewId || view_id2 == kInvalidViewId) {
      continue;
    }
    if (view_id1 > view_id2) {
      SwapCameras(&info);
    }
    view_graph->AddEdge(view_id1, view_id2, info);
  }
}

// Align the reconstructions then evaluate the pose errors.
void EvaluateAlignedPoseError(
    const double position_tolerance_meters,
    const Reconstruction& reference_reconstruction,
    Reconstruction* reconstruction_to_align) {
  // Find the common view names (it should be all views).
  const std::vector<std::string> common_view_names =
      theia::FindCommonViewsByName(reference_reconstruction,
                                   *reconstruction_to_align);
  ASSERT_EQ(common_view_names.size(), reference_reconstruction.NumViews());

  // Align the computed reconstruction with the known ground truth. The ground
  // truth scale is in meters so aligning the computed reconstruction to the
  // ground truth allows us to compute the position error in meters.
  AlignReconstructions(reference_reconstruction, reconstruction_to_align);

  // Ensure the each view's position error is within a tolerance.
  for (int i = 0; i < common_view_names.size(); i++) {
    const ViewId view_id1 =
        reference_reconstruction.ViewIdFromName(common_view_names[i]);
    const ViewId view_id2 =
        reconstruction_to_align->ViewIdFromName(common_view_names[i]);
    const View* view2 = reconstruction_to_align->View(view_id2);
    EXPECT_TRUE(view2->IsEstimated());

    const theia::Camera& camera1 =
        reference_reconstruction.View(view_id1)->Camera();
    const theia::Camera& camera2 = view2->Camera();

    // Compute the position error in meters.
    const double position_error_meters =
        (camera1.GetPosition() - camera2.GetPosition()).norm();
    EXPECT_LT(position_error_meters, position_tolerance_meters);
  }
}

void BuildAndVerifyReconstruction(
    const double position_tolerance_meters,
    const ReconstructionEstimatorOptions& options) {
  ViewGraph view_graph;
  Reconstruction gt_reconstruction, reconstruction;
  ReadInput(&gt_reconstruction, &reconstruction, &view_graph);
  IncrementalReconstructionEstimator reconstruction_estimator(options);

  const ReconstructionEstimatorSummary summary =
      reconstruction_estimator.Estimate(&view_graph, &reconstruction);

  // Ensure all views were estimated.
  EXPECT_EQ(summary.estimated_views.size(), reconstruction.NumViews());

  // Ensure that the reconstruction is somewhat sane.
  EvaluateAlignedPoseError(position_tolerance_meters,
                           gt_reconstruction,
                           &reconstruction);
}

TEST(IncrementalReconstructionEstimator, BasicTest) {
  static const double kPositionToleranceMeters = 1e-2;

  ReconstructionEstimatorOptions options;
  options.reconstruction_estimator_type =
      ReconstructionEstimatorType::INCREMENTAL;
  options.intrinsics_to_optimize = OptimizeIntrinsicsType::NONE;
  BuildAndVerifyReconstruction(kPositionToleranceMeters, options);
}

TEST(IncrementalReconstructionEstimator, RobustCostFunction) {
  static const double kPositionToleranceMeters = 1e-2;

  ReconstructionEstimatorOptions options;
  options.reconstruction_estimator_type =
      ReconstructionEstimatorType::INCREMENTAL;
  options.bundle_adjustment_loss_function_type = LossFunctionType::HUBER;
  options.intrinsics_to_optimize = OptimizeIntrinsicsType::NONE;
  BuildAndVerifyReconstruction(kPositionToleranceMeters, options);
}

TEST(IncrementalReconstructionEstimator, VariableIntrinsics) {
  static const double kPositionToleranceMeters = 1e-2;

  ReconstructionEstimatorOptions options;
  options.reconstruction_estimator_type =
      ReconstructionEstimatorType::INCREMENTAL;
  options.intrinsics_to_optimize = OptimizeIntrinsicsType::FOCAL_LENGTH;
  BuildAndVerifyReconstruction(kPositionToleranceMeters, options);
}

TEST(IncrementalReconstructionEstimator, InitializedReconstruction) {
  static const double kPositionToleranceMeters = 1e-2;

  ReconstructionEstimatorOptions options;
  options.reconstruction_estimator_type =
      ReconstructionEstimatorType::INCREMENTAL;
  options.intrinsics_to_optimize = OptimizeIntrinsicsType::NONE;

  ViewGraph view_graph;
  Reconstruction gt_reconstruction, reconstruction;
  ReadInput(&gt_reconstruction, &reconstruction, &view_graph);

  // Set several of the views to be estimated so that the reconstruction is
  // initialized.
  const std::vector<std::string> initialized_views = {"0000.png", "0001.png",
                                                      "0002.png", "0003.png"};
  for (const std::string& view_name : initialized_views) {
    const ViewId view_id = reconstruction.ViewIdFromName(view_name);
    reconstruction.MutableView(view_id)->SetEstimated(true);
  }
  // Set all tracks as "estimated" that are seen within the initialized views.
  const auto track_ids = reconstruction.TrackIds();
  for (const TrackId track_id : track_ids) {
    reconstruction.MutableTrack(track_id)->SetEstimated(true);
  }
  SetUnderconstrainedTracksToUnestimated(&reconstruction);

  // Estimate the remaining camera and point parameters with incremental SfM.
  IncrementalReconstructionEstimator reconstruction_estimator(options);
  const ReconstructionEstimatorSummary summary =
      reconstruction_estimator.Estimate(&view_graph, &reconstruction);

  // Ensure all views were estimated.
  EXPECT_EQ(summary.estimated_views.size(), reconstruction.NumViews());

  // Ensure that the reconstruction is somewhat sane.
  EvaluateAlignedPoseError(kPositionToleranceMeters,
                           gt_reconstruction,
                           &reconstruction);
}

}  // namespace theia
