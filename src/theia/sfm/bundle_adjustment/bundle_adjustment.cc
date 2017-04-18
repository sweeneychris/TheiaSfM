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

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"

#include <ceres/ceres.h>
#include <glog/logging.h>
#include <algorithm>
#include <memory>
#include <unordered_set>
#include <vector>

#include "theia/util/map_util.h"
#include "theia/util/timer.h"
#include "theia/sfm/bundle_adjustment/create_loss_function.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/create_reprojection_error_cost_function.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/types.h"

namespace theia {

namespace {

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

// Adds camera intrinsic parameters to the problem while optionally holding some
// intrinsics parameters constant.
std::unordered_set<CameraIntrinsicsGroupId> AddCameraIntrinsicsToProblem(
    const OptimizeIntrinsicsType& intrinsics_to_optimize,
    const std::unordered_set<ViewId>& view_ids,
    Reconstruction* reconstruction,
    ceres::Problem* problem) {
  std::unordered_set<CameraIntrinsicsGroupId>
      optimized_camera_intrinsics_groups;
  optimized_camera_intrinsics_groups.reserve(view_ids.size());

  // Loop through all views to optimize and add the intrinsics to the problem.
  for (const ViewId view_id : view_ids) {
    const CameraIntrinsicsGroupId camera_intrinsics_group =
        reconstruction->CameraIntrinsicsGroupIdFromViewId(view_id);

    // If the camera intrinsics group for this camera was already added, skip
    // this view.
    if (ContainsKey(optimized_camera_intrinsics_groups,
                    camera_intrinsics_group)) {
      continue;
    }

    // Otherwise add the camera intrinsics group to the problem.
    optimized_camera_intrinsics_groups.emplace(camera_intrinsics_group);
    std::shared_ptr<CameraIntrinsicsModel>& camera_intrinsics =
        reconstruction->MutableView(view_id)
            ->MutableCamera()
            ->MutableCameraIntrinsics();
    const std::vector<int> constant_intrinsics =
        camera_intrinsics->GetSubsetFromOptimizeIntrinsicsType(
            intrinsics_to_optimize);

    // Set the constant parameters if any are requested.
    if (constant_intrinsics.size() == camera_intrinsics->NumParameters()) {
      problem->AddParameterBlock(camera_intrinsics->mutable_parameters(),
                                 camera_intrinsics->NumParameters());
      problem->SetParameterBlockConstant(
          camera_intrinsics->mutable_parameters());
    } else if (constant_intrinsics.size() > 0) {
      ceres::SubsetParameterization* subset_parameterization =
          new ceres::SubsetParameterization(camera_intrinsics->NumParameters(),
                                            constant_intrinsics);
      problem->AddParameterBlock(camera_intrinsics->mutable_parameters(),
                                 camera_intrinsics->NumParameters(),
                                 subset_parameterization);
    } else {
      problem->AddParameterBlock(camera_intrinsics->mutable_parameters(),
                                 camera_intrinsics->NumParameters());
    }
  }

  return optimized_camera_intrinsics_groups;
}

// Add camera extrinsics to the problem, setting certain values constant based
// on the input settings.
void AddCameraExtrinsicsToProblem(const bool constant_camera_orientation,
                                  const bool constant_camera_position,
                                  const std::unordered_set<ViewId>& view_ids,
                                  Reconstruction* reconstruction,
                                  ceres::Problem* problem) {
  if (constant_camera_orientation && constant_camera_position) {
    for (const ViewId view_id : view_ids) {
      View* view = reconstruction->MutableView(view_id);
      Camera* camera = view->MutableCamera();
      problem->AddParameterBlock(camera->mutable_extrinsics(),
                                 Camera::kExtrinsicsSize);
      problem->SetParameterBlockConstant(camera->mutable_extrinsics());

    }
  } else if (!constant_camera_orientation && !constant_camera_position) {
    for (const ViewId view_id : view_ids) {
      View* view = reconstruction->MutableView(view_id);
      Camera* camera = view->MutableCamera();
      problem->AddParameterBlock(camera->mutable_extrinsics(),
                                 Camera::kExtrinsicsSize);
    }
  } else {
    std::vector<int> constant_extrinsics;
    if (constant_camera_position) {
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
    for (const ViewId view_id : view_ids) {
      View* view = reconstruction->MutableView(view_id);
      Camera* camera = view->MutableCamera();
      problem->AddParameterBlock(camera->mutable_extrinsics(),
                                 Camera::kExtrinsicsSize,
                                 subset_parameterization);
    }
  }
}

}  // namespace

// Bundle adjust the entire model.
BundleAdjustmentSummary BundleAdjustPartialReconstruction(
    const BundleAdjustmentOptions& options,
    const std::unordered_set<ViewId>& view_ids,
    const std::unordered_set<TrackId>& track_ids,
    Reconstruction* reconstruction) {
  CHECK_NOTNULL(reconstruction);
  BundleAdjustmentSummary summary;
  static const int kTrackSize = 4;

  // Start setup timer.
  Timer timer;

  // Get the loss function that will be used for BA.
  ceres::Problem::Options problem_options;
  std::unique_ptr<ceres::LossFunction> loss_function =
      CreateLossFunction(options.loss_function_type, options.robust_loss_width);
  problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;;
  ceres::Problem problem(problem_options);

  // Set solver options.
  ceres::Solver::Options solver_options;
  SetSolverOptions(options, &solver_options);
  ceres::ParameterBlockOrdering* parameter_ordering =
      solver_options.linear_solver_ordering.get();

  // Add all camera intrinsics to the problem. We keep track of which camera
  // intrinsics groups are optimized so that later on we can determine which
  // intrinsics need to be constant or variable during the optimization.
  const std::unordered_set<CameraIntrinsicsGroupId>
      optimized_camera_intrinsics_groups = AddCameraIntrinsicsToProblem(
          options.intrinsics_to_optimize, view_ids, reconstruction, &problem);

  // Add all of the camera extrinsics to the problem, setting orientation and/or
  // position constant as appropriate.
  AddCameraExtrinsicsToProblem(options.constant_camera_orientation,
                               options.constant_camera_position,
                               view_ids,
                               reconstruction,
                               &problem);

  // Per recommendation of Ceres documentation we group the parameters by points
  // (group 0) and camera parameters (group 1) so that the points are eliminated
  // first then the cameras.
  for (const ViewId view_id : view_ids) {
    View* view = CHECK_NOTNULL(reconstruction->MutableView(view_id));
    // Only optimize estimated views.
    if (!view->IsEstimated()) {
      continue;
    }

    // Fetch the camera that will be optimized.
    Camera* camera = view->MutableCamera();

    // Get a pointer to the shared camera intrinsics.
    std::shared_ptr<CameraIntrinsicsModel>& shared_intrinsics =
        camera->MutableCameraIntrinsics();

    // Add camera parameters to groups 1 and 2. The extrinsics *must* belong to
    // group 2. This is because inner iterations uses a reverse ordering of
    // elimination and the Schur-based solvers require the first group to be an
    // independent set. Since the intrinsics may be shared, they are not
    // guaranteed to form an independent set and so we must use the extrinsics
    // in group 2.
    parameter_ordering->AddElementToGroup(camera->mutable_extrinsics(), 2);
    parameter_ordering->AddElementToGroup(camera->mutable_intrinsics(), 1);

    // Add residuals for all tracks in the view.
    for (const TrackId track_id : view->TrackIds()) {
      const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
      Track* track = CHECK_NOTNULL(reconstruction->MutableTrack(track_id));
      // Only consider tracks with an estimated 3d point.
      if (!track->IsEstimated()) {
        continue;
      }

      problem.AddResidualBlock(
          CreateReprojectionErrorCostFunction(
              camera->GetCameraIntrinsicsModelType(), *feature),
          loss_function.get(),
          camera->mutable_extrinsics(),
          camera->mutable_intrinsics(),
          track->MutablePoint()->data());
      // Add the point to group 0.
      parameter_ordering->AddElementToGroup(track->MutablePoint()->data(), 0);
      problem.SetParameterBlockConstant(track->MutablePoint()->data());
    }
  }

  // The previous loop gives us residuals for all tracks in all the views that
  // we want to optimize. However, the tracks should still be constrained by
  // *all* views that observe it, not just the ones we want to optimize. Here,
  // we add in any views that were not part of the first loop and we keep them
  // constant during the optimization. We need to take care to ensure that
  std::unordered_set<ViewId> potentially_constant_camera_intrinsics_groups;
  for (const TrackId track_id : track_ids) {
    Track* track = CHECK_NOTNULL(reconstruction->MutableTrack(track_id));
    if (!track->IsEstimated()) {
      continue;
    }

    // Add the track to the problem and set the parameter ordering for Schur
    // elimination.
    problem.AddParameterBlock(track->MutablePoint()->data(), kTrackSize);
    parameter_ordering->AddElementToGroup(track->MutablePoint()->data(), 0);
    problem.SetParameterBlockVariable(track->MutablePoint()->data());

    const auto& observed_view_ids = track->ViewIds();
    for (const ViewId view_id : observed_view_ids) {
      View* view = CHECK_NOTNULL(reconstruction->MutableView(view_id));
      // Only optimize estimated views that have not already been added.
      if (ContainsKey(view_ids, view_id) || !view->IsEstimated()) {
        continue;
      }

      const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
      Camera* camera = view->MutableCamera();

      // Add the residual for the track to the problem. The shared intrinsics
      // parameter block will be set to constant after the loop if no optimized
      // cameras share the same camera intrinsics.
      problem.AddResidualBlock(
          CreateReprojectionErrorCostFunction(
              camera->GetCameraIntrinsicsModelType(), *feature),
          loss_function.get(),
          camera->mutable_extrinsics(),
          camera->mutable_intrinsics(),
          track->MutablePoint()->data());

      // Add camera parameters to groups 1 and 2.
      parameter_ordering->AddElementToGroup(camera->mutable_extrinsics(), 2);
      parameter_ordering->AddElementToGroup(camera->mutable_intrinsics(), 1);

      // Any camera that reaches this point was not part of the first loop, so
      // we do not want to optimize it.
      problem.SetParameterBlockConstant(camera->mutable_extrinsics());

      // Mark the camera intrinsics as "potentially constant." We only set the
      // parameter block to constant if the shared intrinsics are not shared
      // with cameras that are being optimized.
      const CameraIntrinsicsGroupId intrinsics_group_id =
        reconstruction->CameraIntrinsicsGroupIdFromViewId(view_id);
      potentially_constant_camera_intrinsics_groups.emplace(
          intrinsics_group_id);
    }
  }

  // Set any potentially constant camera intrinsics groups that do not have a
  // mutable view to be constant.
  for (const CameraIntrinsicsGroupId shared_intrinsics_id :
       potentially_constant_camera_intrinsics_groups) {
    // If the camera intrinsics group does not contain a view that is being
    // optimized, then we set the intrinsics to be constant.
    if (!ContainsKey(optimized_camera_intrinsics_groups,
                     shared_intrinsics_id)) {
      // First, fetch any view from the camera intrinsics group.
      const ViewId view_id_in_intrinsics_group =
          *reconstruction->GetViewsInCameraIntrinsicGroup(shared_intrinsics_id)
               .begin();
      // Then, set the intrinsics to be constant. Since the intrinsics share the
      // same memory for the whole group, this will set the entire camera
      // intrinsics group to be constant.
      Camera* intrinsics_group_camera =
          reconstruction->MutableView(view_id_in_intrinsics_group)
              ->MutableCamera();
      problem.SetParameterBlockConstant(
          intrinsics_group_camera->mutable_intrinsics());
    }
  }

  // NOTE: cmsweeney found a thread on the Ceres Solver email group that
  // indicated using the reverse BA order (i.e., using cameras then points) is a
  // good idea for inner iterations.
  if (solver_options.use_inner_iterations) {
    solver_options.inner_iteration_ordering.reset(
        new ceres::ParameterBlockOrdering(*parameter_ordering));
    solver_options.inner_iteration_ordering->Reverse();
  }

  // Solve the problem.
  const double internal_setup_time = timer.ElapsedTimeInSeconds();
  ceres::Solver::Summary solver_summary;
  ceres::Solve(solver_options, &problem, &solver_summary);
  LOG_IF(INFO, options.verbose) << solver_summary.FullReport();

  // Set the BundleAdjustmentSummary.
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

// Bundle adjust the specified views and all tracks observed by those views.
BundleAdjustmentSummary BundleAdjustReconstruction(
    const BundleAdjustmentOptions& options,
    Reconstruction* reconstruction) {
  const auto& view_ids = reconstruction->ViewIds();
  const auto& track_ids = reconstruction->TrackIds();
  const std::unordered_set<ViewId> view_ids_set(view_ids.begin(),
                                                view_ids.end());
  const std::unordered_set<TrackId> track_ids_set(track_ids.begin(),
                                                  track_ids.end());
  return BundleAdjustPartialReconstruction(options,
                                           view_ids_set,
                                           track_ids_set,
                                           reconstruction);
}

// Bundle adjust a single view.
BundleAdjustmentSummary BundleAdjustView(const BundleAdjustmentOptions& options,
                                         const ViewId view_id,
                                         Reconstruction* reconstruction) {
  std::unordered_set<ViewId> view_ids;
  view_ids.insert(view_id);
  std::unordered_set<TrackId> track_ids;
  BundleAdjustmentOptions ba_options = options;
  ba_options.linear_solver_type = ceres::DENSE_QR;
  ba_options.use_inner_iterations = false;
  return BundleAdjustPartialReconstruction(ba_options,
                                           view_ids,
                                           track_ids,
                                           reconstruction);
}

// Bundle adjust a single track.
BundleAdjustmentSummary BundleAdjustTrack(
    const BundleAdjustmentOptions& options,
    const TrackId my_track_id,
    Reconstruction* reconstruction) {
  std::unordered_set<ViewId> view_ids;
  std::unordered_set<TrackId> track_ids;
  track_ids.insert(my_track_id);
  BundleAdjustmentOptions ba_options = options;
  ba_options.linear_solver_type = ceres::DENSE_QR;
  ba_options.use_inner_iterations = false;
  return BundleAdjustPartialReconstruction(ba_options,
                                           view_ids,
                                           track_ids,
                                           reconstruction);
}

}  // namespace theia
