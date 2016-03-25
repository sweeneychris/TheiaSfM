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
#include <memory>
#include <unordered_set>
#include <vector>

#include "theia/util/map_util.h"
#include "theia/util/timer.h"
#include "theia/sfm/bundle_adjustment/create_loss_function.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/reprojection_error.h"
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

// Determine which camera intrinsics to optimize by investigating the individual
// bits of intrinsics_to_optimize.
std::vector<int> GetIntrinsicsToOptimize(
    const OptimizeIntrinsicsType& intrinsics_to_optimize) {
  std::vector<int> constant_intrinsics;
  if (intrinsics_to_optimize == OptimizeIntrinsicsType::ALL) {
    return constant_intrinsics;
  }

  if ((intrinsics_to_optimize &
      OptimizeIntrinsicsType::FOCAL_LENGTH) == OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(Camera::kExtrinsicsSize +
                                     Camera::FOCAL_LENGTH);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::ASPECT_RATIO) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(Camera::kExtrinsicsSize +
                                     Camera::ASPECT_RATIO);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::SKEW) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(Camera::kExtrinsicsSize + Camera::SKEW);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::PRINCIPAL_POINTS) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(Camera::kExtrinsicsSize +
                                     Camera::PRINCIPAL_POINT_X);
    constant_intrinsics.emplace_back(Camera::kExtrinsicsSize +
                                     Camera::PRINCIPAL_POINT_Y);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::RADIAL_DISTORTION) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(Camera::kExtrinsicsSize +
                                     Camera::RADIAL_DISTORTION_1);
    constant_intrinsics.emplace_back(Camera::kExtrinsicsSize +
                                     Camera::RADIAL_DISTORTION_2);
  }
  return constant_intrinsics;
}

// Adds camera parameters to the problem while optionally holding some
// intrinsics parameters constant.
void AddCameraParametersToProblem(const std::vector<int>& constant_intrinsics,
                                  double* camera_parameters,
                                  ceres::Problem* problem) {
  if (constant_intrinsics.size() > 0) {
    ceres::SubsetParameterization* subset_parameterization =
      new ceres::SubsetParameterization(Camera::kParameterSize,
                                        constant_intrinsics);
    problem->AddParameterBlock(camera_parameters,
                               Camera::kParameterSize,
                               subset_parameterization);
  } else {
    problem->AddParameterBlock(camera_parameters, Camera::kParameterSize);
  }

  // Set bounds for certain camera parameters to make sure they are reasonable.
  problem->SetParameterLowerBound(
      camera_parameters, Camera::kExtrinsicsSize + Camera::FOCAL_LENGTH, 0.0);
  problem->SetParameterLowerBound(
      camera_parameters, Camera::kExtrinsicsSize + Camera::ASPECT_RATIO, 0.0);
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

  // Obtain which params will be constant during optimization.
  const std::vector<int> constant_intrinsics =
      GetIntrinsicsToOptimize(options.intrinsics_to_optimize);

  // Per recommendation of Ceres documentation we group the parameters by points
  // (group 0) and camera parameters (group 1) so that the points are eliminated
  // first then the cameras.
  for (const ViewId view_id : view_ids) {
    View* view = CHECK_NOTNULL(reconstruction->MutableView(view_id));
    // Only optimize estimated views.
    if (!view->IsEstimated()) {
      continue;
    }

    Camera* camera = view->MutableCamera();
    // This function will add all camera parameters to the problem and will keep
    // the intrinsic params constant if desired.
    AddCameraParametersToProblem(constant_intrinsics,
                                 camera->mutable_parameters(),
                                 &problem);

    // Add camera parameters to group 1.
    parameter_ordering->AddElementToGroup(camera->mutable_parameters(), 1);

    // Add residuals for all tracks in the view.
    for (const TrackId track_id : view->TrackIds()) {
      const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
      Track* track = CHECK_NOTNULL(reconstruction->MutableTrack(track_id));
      // Only consider tracks with an estimated 3d point.
      if (!track->IsEstimated()) {
        continue;
      }

      problem.AddResidualBlock(
          ReprojectionError::Create(*feature),
          loss_function.get(),
          camera->mutable_parameters(),
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
  // constant during the optimization.
  for (const TrackId track_id : track_ids) {
    Track* track = CHECK_NOTNULL(reconstruction->MutableTrack(track_id));
    if (!track->IsEstimated()) {
      continue;
    }

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

      Camera* camera = view->MutableCamera();
      const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
      problem.AddResidualBlock(
          ReprojectionError::Create(*feature),
          loss_function.get(),
          camera->mutable_parameters(),
          track->MutablePoint()->data());

      // Add camera parameters to group 1.
      parameter_ordering->AddElementToGroup(camera->mutable_parameters(), 1);
      // Any camera that reaches this point was not part of the first loop, so
      // we do not want to optimize it.
      problem.SetParameterBlockConstant(camera->mutable_parameters());
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
                                         Reconstruction* reconstruction){
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
BundleAdjustmentSummary BundleAdjustTrack(const BundleAdjustmentOptions& options,
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
