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

#include "theia/sfm/bundle_adjustment/bundle_adjuster.h"

#include <algorithm>
#include <ceres/ceres.h>
#include <glog/logging.h>
#include <memory>
#include <unordered_set>
#include <vector>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/bundle_adjustment/create_loss_function.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/create_reprojection_error_cost_function.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/types.h"
#include "theia/util/map_util.h"
#include "theia/util/timer.h"

namespace theia {

namespace {

// Set the solver options to defaults.
void SetSolverOptions(const BundleAdjustmentOptions& options,
                      ceres::Solver::Options* solver_options) {
  solver_options->linear_solver_type = options.linear_solver_type;
  solver_options->preconditioner_type = options.preconditioner_type;
  solver_options->visibility_clustering_type =
      options.visibility_clustering_type;
  solver_options->logging_type =
      options.verbose ? ceres::PER_MINIMIZER_ITERATION : ceres::SILENT;
  solver_options->num_threads = options.num_threads;
  solver_options->num_linear_solver_threads = options.num_threads;
  solver_options->max_num_iterations = options.max_num_iterations;
  solver_options->max_solver_time_in_seconds =
      options.max_solver_time_in_seconds;
  solver_options->use_inner_iterations = options.use_inner_iterations;
  solver_options->function_tolerance = options.function_tolerance;
  solver_options->gradient_tolerance = options.gradient_tolerance;
  solver_options->parameter_tolerance = options.parameter_tolerance;
  solver_options->max_trust_region_radius = options.max_trust_region_radius;

  // Solver options takes ownership of the ordering so that we can order the BA
  // problem by points and cameras.
  solver_options->linear_solver_ordering.reset(
      new ceres::ParameterBlockOrdering);
}

}  // namespace

BundleAdjuster::BundleAdjuster(const BundleAdjustmentOptions& options,
                               Reconstruction* reconstruction)
    : options_(options), reconstruction_(reconstruction) {
  CHECK_NOTNULL(reconstruction);

  // Start setup timer.
  timer_.Reset();

  // Get the loss function that will be used for BA.
  loss_function_ =
      CreateLossFunction(options.loss_function_type, options.robust_loss_width);
  ceres::Problem::Options problem_options;
  problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
  problem_.reset(new ceres::Problem(problem_options));

  // Set solver options.
  SetSolverOptions(options, &solver_options_);
  parameter_ordering_ = solver_options_.linear_solver_ordering.get();
}

void BundleAdjuster::AddView(const ViewId view_id) {
  View* view = CHECK_NOTNULL(reconstruction_->MutableView(view_id));

  // Only optimize estimated views.
  if (!view->IsEstimated() || ContainsKey(optimized_views_, view_id)) {
    return;
  }

  // Mark the view as optimized.
  optimized_views_.emplace(view_id);

  // Mark the camera intrinsics as optimized.
  const CameraIntrinsicsGroupId intrinsics_group_id =
      reconstruction_->CameraIntrinsicsGroupIdFromViewId(view_id);
  optimized_camera_intrinsics_groups_.emplace(intrinsics_group_id);

  // Fetch the camera that will be optimized.
  Camera* camera = view->MutableCamera();

  // Add camera parameters to groups 1 and 2. The extrinsics *must* belong to
  // group 2. This is because inner iterations uses a reverse ordering of
  // elimination and the Schur-based solvers require the first group to be an
  // independent set. Since the intrinsics may be shared, they are not
  // guaranteed to form an independent set and so we must use the extrinsics
  // in group 2.
  parameter_ordering_->AddElementToGroup(camera->mutable_extrinsics(),
                                         kExtrinsicsParameterGroup);
  parameter_ordering_->AddElementToGroup(camera->mutable_intrinsics(),
                                         kIntrinsicsParameterGroup);

  // Add residuals for all tracks in the view.
  for (const TrackId track_id : view->TrackIds()) {
    const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
    Track* track = CHECK_NOTNULL(reconstruction_->MutableTrack(track_id));
    // Only consider tracks with an estimated 3d point.
    if (!track->IsEstimated()) {
      continue;
    }

    // Add the reprojection error to the optimization.
    AddReprojectionErrorResidual(*feature, camera, track);

    // Add the point to group 0.
    parameter_ordering_->AddElementToGroup(track->MutablePoint()->data(),
                                           kTrackParameterGroup);
    problem_->SetParameterBlockConstant(track->MutablePoint()->data());
  }
}

void BundleAdjuster::AddTrack(const TrackId track_id) {
  Track* track = CHECK_NOTNULL(reconstruction_->MutableTrack(track_id));
  // Only optimize estimated tracks.
  if (!track->IsEstimated() || ContainsKey(optimized_tracks_, track_id)) {
    return;
  }

  // Mark the track as optimized.
  optimized_tracks_.emplace(track_id);

  // Add all observations of the track to the problem.
  const auto& observed_view_ids = track->ViewIds();
  for (const ViewId view_id : observed_view_ids) {
    View* view = CHECK_NOTNULL(reconstruction_->MutableView(view_id));
    // Only optimize estimated views that have not already been added.
    if (ContainsKey(optimized_views_, view_id) || !view->IsEstimated()) {
      continue;
    }

    const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
    Camera* camera = view->MutableCamera();

    // Add the reprojection error to the optimization.
    AddReprojectionErrorResidual(*feature, camera, track);

    // Add camera parameters to groups for Schur elimination.
    parameter_ordering_->AddElementToGroup(camera->mutable_extrinsics(),
                                           kExtrinsicsParameterGroup);
    parameter_ordering_->AddElementToGroup(camera->mutable_intrinsics(),
                                           kIntrinsicsParameterGroup);

    // Any camera that reaches this point was not added by AddView() and so we
    // want to mark it as constant.
    problem_->SetParameterBlockConstant(camera->mutable_extrinsics());

    // Mark the camera intrinsics as "potentially constant." We only set the
    // parameter block to constant if the shared intrinsics are not shared
    // with cameras that are being optimized.
    const CameraIntrinsicsGroupId intrinsics_group_id =
        reconstruction_->CameraIntrinsicsGroupIdFromViewId(view_id);
    potentially_constant_camera_intrinsics_groups_.emplace(intrinsics_group_id);
  }

  // Set the parameter ordering for Schur elimination. We do this after the loop
  // above so that the track is already added to the problem.
  parameter_ordering_->AddElementToGroup(track->MutablePoint()->data(),
                                         kTrackParameterGroup);
  problem_->SetParameterBlockVariable(track->MutablePoint()->data());
}

BundleAdjustmentSummary BundleAdjuster::Optimize() {
  // Set extrinsics parameterization of the camera poses. This will set
  // orientation and/or positions as constant if desired.
  SetCameraExtrinsicsParameterization();

  // Set intrinsics group parameterization. This will control which of the
  // intrinsics parameters or optimized or held constant. Note that each camera
  // intrinsics group may be a different camera and/or a different camera
  // intrinsics model.
  SetCameraIntrinsicsParameterization();

  // NOTE: csweeney found a thread on the Ceres Solver email group that
  // indicated using the reverse BA order (i.e., using cameras then points) is a
  // good idea for inner iterations.
  if (solver_options_.use_inner_iterations) {
    solver_options_.inner_iteration_ordering.reset(
        new ceres::ParameterBlockOrdering(*parameter_ordering_));
    solver_options_.inner_iteration_ordering->Reverse();
  }

  // Solve the problem.
  const double internal_setup_time = timer_.ElapsedTimeInSeconds();
  ceres::Solver::Summary solver_summary;
  ceres::Solve(solver_options_, problem_.get(), &solver_summary);
  LOG_IF(INFO, options_.verbose) << solver_summary.FullReport();

  // Set the BundleAdjustmentSummary.
  BundleAdjustmentSummary summary;
  summary.setup_time_in_seconds =
      internal_setup_time + solver_summary.preprocessor_time_in_seconds;
  summary.solve_time_in_seconds = solver_summary.total_time_in_seconds;
  summary.initial_cost = solver_summary.initial_cost;
  summary.final_cost = solver_summary.final_cost;

  // This only indicates whether the optimization was successfully run and makes
  // no guarantees on the quality or convergence.
  summary.success = solver_summary.IsSolutionUsable();

  return summary;
}

void BundleAdjuster::SetCameraExtrinsicsParameterization() {
  // If all camera extrinsics are mutable then simply return since all
  // parameters have already been added to the problem as mutable.
  if (!options_.constant_camera_orientation &&
      !options_.constant_camera_position) {
    return;
  } else if (options_.constant_camera_orientation &&
             options_.constant_camera_position) {
    // If all extrinsics are constant then mark the entire parameter block
    // as constant.
    for (const ViewId view_id : optimized_views_) {
      View* view = reconstruction_->MutableView(view_id);
      Camera* camera = view->MutableCamera();
      problem_->SetParameterBlockConstant(camera->mutable_extrinsics());
    }
  } else {
    // Otherwise, only mark the relevant parameters as constant.
    std::vector<int> constant_extrinsics;
    if (options_.constant_camera_position) {
      constant_extrinsics.emplace_back(Camera::POSITION + 0);
      constant_extrinsics.emplace_back(Camera::POSITION + 1);
      constant_extrinsics.emplace_back(Camera::POSITION + 2);
    } else {
      constant_extrinsics.emplace_back(Camera::ORIENTATION + 0);
      constant_extrinsics.emplace_back(Camera::ORIENTATION + 1);
      constant_extrinsics.emplace_back(Camera::ORIENTATION + 2);
    }
    ceres::SubsetParameterization* subset_parameterization =
        new ceres::SubsetParameterization(Camera::kExtrinsicsSize,
                                          constant_extrinsics);
    for (const ViewId view_id : optimized_views_) {
      View* view = reconstruction_->MutableView(view_id);
      Camera* camera = view->MutableCamera();
      problem_->SetParameterization(camera->mutable_extrinsics(),
                                    subset_parameterization);
    }
  }
}

void BundleAdjuster::SetCameraIntrinsicsParameterization() {
  // Loop through all optimized camera intrinsics groups to set the intrinsics
  // parameterization.
  for (const CameraIntrinsicsGroupId camera_intrinsics_group :
       optimized_camera_intrinsics_groups_) {
    // Get a representative view for the intrinsics group.
    const auto& views_in_intrinsics_groups =
        reconstruction_->GetViewsInCameraIntrinsicGroup(
            camera_intrinsics_group);
    CHECK(!views_in_intrinsics_groups.empty());
    const ViewId representative_view_id = *views_in_intrinsics_groups.begin();

    // Get the shared camera intrinsics parameters.
    std::shared_ptr<CameraIntrinsicsModel>& camera_intrinsics =
        reconstruction_->MutableView(representative_view_id)
            ->MutableCamera()
            ->MutableCameraIntrinsics();

    // Get the subset parameterization of the intrinsics to keep constant.
    const std::vector<int> constant_intrinsics =
        camera_intrinsics->GetSubsetFromOptimizeIntrinsicsType(
            options_.intrinsics_to_optimize);

    // Set the constant parameters if any are requested.
    if (constant_intrinsics.size() == camera_intrinsics->NumParameters()) {
      problem_->SetParameterBlockConstant(
          camera_intrinsics->mutable_parameters());
    } else if (constant_intrinsics.size() > 0) {
      ceres::SubsetParameterization* subset_parameterization =
          new ceres::SubsetParameterization(camera_intrinsics->NumParameters(),
                                            constant_intrinsics);
      problem_->SetParameterization(camera_intrinsics->mutable_parameters(),
                                    subset_parameterization);
    }
  }

  // Set camera intrinsics to be constant if no cameras in the intrinsics group
  // are being optimized.
  for (const CameraIntrinsicsGroupId camera_intrinsics_group :
       potentially_constant_camera_intrinsics_groups_) {
    // If the intrinsics group is being optimized, do not mark it as constant.
    if (ContainsKey(optimized_camera_intrinsics_groups_,
                    camera_intrinsics_group)) {
      continue;
    }

    // Get a representative view for the intrinsics group.
    const auto& views_in_intrinsics_groups =
        reconstruction_->GetViewsInCameraIntrinsicGroup(
            camera_intrinsics_group);
    CHECK(!views_in_intrinsics_groups.empty());
    const ViewId representative_view_id = *views_in_intrinsics_groups.begin();

    // Get the shared camera intrinsics parameters.
    std::shared_ptr<CameraIntrinsicsModel>& camera_intrinsics =
        reconstruction_->MutableView(representative_view_id)
            ->MutableCamera()
            ->MutableCameraIntrinsics();

    // Set the intrinsics to be constant.
    problem_->SetParameterBlockConstant(
        camera_intrinsics->mutable_parameters());
  }
}

void BundleAdjuster::AddReprojectionErrorResidual(const Feature& feature,
                                                  Camera* camera,
                                                  Track* track) {
  // Add the residual for the track to the problem. The shared intrinsics
  // parameter block will be set to constant after the loop if no optimized
  // cameras share the same camera intrinsics.
  problem_->AddResidualBlock(
      CreateReprojectionErrorCostFunction(
          camera->GetCameraIntrinsicsModelType(), feature),
      loss_function_.get(),
      camera->mutable_extrinsics(),
      camera->mutable_intrinsics(),
      track->MutablePoint()->data());
}

}  // namespace theia
