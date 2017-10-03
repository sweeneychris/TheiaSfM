// Copyright (C) 2017 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#ifndef THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUSTER_H_
#define THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUSTER_H_

#include <ceres/ceres.h>
#include <ceres/types.h>
#include <unordered_set>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/types.h"
#include "theia/util/timer.h"

namespace theia {
class Camera;
class Reconstruction;
class Track;

// This class sets up nonlinear optimization problems for bundle adjustment.
// Bundle adjustment problems are set up by adding views and tracks to be
// optimized. Only the views and tracks supplied with AddView and AddTrack will
// be optimized. All other parameters are held constant.
//
// NOTE: It is required that AddViews is called before AddTracks if any views
// are being optimized.
class BundleAdjuster {
 public:
  // Set up the bundle adjuster. The reconstruction will be modified during
  // bundle adjustment.
  BundleAdjuster(const BundleAdjustmentOptions& options,
                 Reconstruction* reconstruction);

  // Add a view to be optimized with bundle adjustment. A residual is created
  // for each estimated track that the view observes.
  virtual void AddView(const ViewId view_id);

  // Add a track to be optimized with bundle adjustment. A residual is created
  // for each estimated view that observes the track.
  virtual void AddTrack(const TrackId track_id);

  // After AddView and AddTrack have been called, optimize the provided views
  // and tracks with bundle adjustment.
  virtual BundleAdjustmentSummary Optimize();

 private:
  static const int kTrackParameterGroup = 0;
  static const int kIntrinsicsParameterGroup = 1;
  static const int kExtrinsicsParameterGroup = 2;

  // Add all camera extrinsics and intrinsics to the optimization problem.
  void SetCameraExtrinsicsParameterization();
  void SetCameraIntrinsicsParameterization();

  // Add the reprojection error residual to the problem.
  virtual void AddReprojectionErrorResidual(const Feature& feature,
                                            Camera* camera,
                                            Track* track);

  const BundleAdjustmentOptions options_;
  Reconstruction* reconstruction_;
  Timer timer_;

  // Ceres problem for optimization.
  std::unique_ptr<ceres::Problem> problem_;
  ceres::Solver::Options solver_options_;

  // The potentially robust loss function to use for reprojection error
  // minimization.
  std::unique_ptr<ceres::LossFunction> loss_function_;
  // The parameter group ordering for bundle adjustment.
  ceres::ParameterBlockOrdering* parameter_ordering_;

  // The optimized views.
  std::unordered_set<ViewId> optimized_views_;
  // The optimized tracks.
  std::unordered_set<TrackId> optimized_tracks_;

  // The intrinsics groups that are optimized.
  std::unordered_set<CameraIntrinsicsGroupId>
      optimized_camera_intrinsics_groups_;

  // Intrinsics groups that have at least 1 camera marked as "const" during
  // optimization. Only the intrinsics that have no optimized cameras are kept
  // as constant during optimization.
  std::unordered_set<CameraIntrinsicsGroupId>
      potentially_constant_camera_intrinsics_groups_;
};

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUSTER_H_
