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

#ifndef THEIA_SFM_RECONSTRUCTION_ESTIMATOR_OPTIONS_H_
#define THEIA_SFM_RECONSTRUCTION_ESTIMATOR_OPTIONS_H_

namespace theia {

enum class ReconstructionEstimatorType {
  NONLINEAR = 0,
};

// Options for the reconstruction estimation.
struct ReconstructionEstimatorOptions {
  // Type of reconstruction estimation to use.
  ReconstructionEstimatorType reconstruction_estimator_type =
      ReconstructionEstimatorType::NONLINEAR;

  // Number of threads to use.
  int num_threads = 1;

  // Maximum reprojection error. This is used to determine inlier
  // correspondences for absolute pose estimation. Additionally, this is the
  // threshold used for filtering outliers after bundle adjustment.
  double max_reprojection_error_in_pixels = 5.0;

  // Any edges in the view graph with fewer than min_num_two_view_inliers will
  // be removed as an initial filtering step.
  int min_num_two_view_inliers = 30;

  // After computing a model and performing an initial BA, the reconstruction
  // can be further improved (and even densified) if we attempt (again) to
  // retriangulate any tracks that are currently unestimated. For each
  // retriangulation iteration we do the following:
  //   1. Remove features that are above max_reprojection_error_in_pixels.
  //   2. Triangulate all unestimated tracks.
  //   3. Perform full bundle adjustment.
  int num_retriangulation_iterations = 1;

  // --------------- RANSAC Options --------------- //
  double ransac_confidence = 0.9999;
  int ransac_min_iterations = 50;
  int ransac_max_iterations = 1000;
  bool ransac_use_mle = true;

  // --------------- Rotation Filtering Options --------------- //

  // After orientations are estimated, view pairs may be filtered/removed if the
  // relative rotation of the view pair differs from the relative rotation
  // formed by the global orientation estimations. Adjust this threshold to
  // control the threshold at which rotations are filtered. See
  // theia/sfm/filter_view_pairs_from_orientation.h
  double rotation_filtering_max_difference_degrees = 5.0;

  // --------------- Position Filtering Options --------------- //

  // Refine the relative translations based on the epipolar error and known
  // rotation estimations. This improve the quality of the translation
  // estimation.
  bool refine_relative_translations_after_rotation_estimation = true;

  // Before the camera positions are estimated, it is wise to remove any
  // relative translations estimates that are low quality. See
  // theia/sfm/filter_view_pairs_from_relative_translation.h
  int translation_filtering_num_iterations = 48;
  double translation_filtering_projection_tolerance = 0.1;

  // --------------- Nonlinear Rotation Estimation Options --------------- //

  // Robust loss function scales for nonlinear estimation.
  double rotation_estimation_robust_loss_scale = 0.1;

  // --------------- Nonlinear Position Estimation Options --------------- //
  double position_estimation_robust_loss_scale = 1.0;

  // Number of point to camera correspondences used for nonlinear position
  // estimation.
  int position_estimation_min_num_tracks_per_view = 10;
  // Weight of point to camera constraints with respect to camera to camera
  // constraints.
  double position_estimation_point_to_camera_weight = 0.5;

  // --------------- Triangulation Options --------------- //

  // Minimum angle required between a 3D point and 2 viewing rays in order to
  // consider triangulation a success.
  double min_triangulation_angle_degrees = 3.0;

  // The reprojection error to use for determining valid triangulation.
  double triangulation_max_reprojection_error_in_pixels = 10.0;
  // Bundle adjust a track immediately after estimating it.
  bool bundle_adjust_tracks = true;

  // --------------- Bundle Adjustment Options --------------- //

  // Use SPARSE_SCHUR for problems smaller than this size and ITERATIVE_SCHUR
  // for problems larger than this size.
  int min_cameras_for_iterative_solver = 1000;

  // If accurate calibration is known ahead of time then it is recommended to
  // set the camera intrinsics constant during bundle adjustment.
  bool constant_camera_intrinsics = false;
};

}  // namespace theia

#endif  // THEIA_SFM_RECONSTRUCTION_ESTIMATOR_OPTIONS_H_
