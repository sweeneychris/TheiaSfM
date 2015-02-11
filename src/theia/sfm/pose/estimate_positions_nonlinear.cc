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

#include "theia/sfm/pose/estimate_positions_nonlinear.h"

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Core>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "theia/util/map_util.h"
#include "theia/util/util.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/pose/pairwise_translation_error.h"
#include "theia/sfm/types.h"

namespace theia {
namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

Vector3d GetRotatedTranslation(const Vector3d& rotation_angle_axis,
                               const Vector3d& translation) {
  Matrix3d rotation;
  ceres::AngleAxisToRotationMatrix(
      rotation_angle_axis.data(),
      ceres::ColumnMajorAdapter3x3(rotation.data()));
  return rotation.transpose() * translation;
}

Vector3d GetRotatedFeatureRay(const Camera& camera,
                              const Vector3d& orientation,
                              const Feature& feature) {
  Camera temp_camera = camera;
  temp_camera.SetOrientationFromAngleAxis(orientation);
  // Get the image ray rotated into the world reference frame.
  return camera.PixelToUnitDepthRay(feature).normalized();
}

// Sorts the pairs such that the number of views (i.e. the int) is sorted in
// descending order.
bool CompareViewsPerTrack(const std::pair<TrackId, int>& t1,
                          const std::pair<TrackId, int>& t2) {
  return t1.second > t2.second;
}

}  // namespace

NonlinearPositionEstimator::NonlinearPositionEstimator(
    const NonlinearPositionEstimatorOptions& options,
    const Reconstruction& reconstruction,
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs)
    : options_(options),
      reconstruction_(reconstruction),
      view_pairs_(view_pairs) {
  CHECK_GT(options_.num_threads, 0);
  CHECK_GE(options_.min_num_points_per_view, 0);
  CHECK_GT(options_.point_to_camera_weight, 0);
  CHECK_GT(options_.robust_loss_width, 0);
  CHECK_GE(options_.max_num_reweighted_iterations, 0);
  CHECK_GT(options_.reweighted_convergence_tolerance, 0);
}

bool NonlinearPositionEstimator::EstimatePositions(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  CHECK_NOTNULL(positions);
  triangulated_points_.clear();
  problem_.reset(new ceres::Problem());

  // Iterative schur is only used if the problem is large enough, otherwise
  // sparse schur is used.
  static const int kMinNumCamerasForIterativeSchur = 1000;

  // Initialize positions to be random.
  InitializeRandomPositions(orientations, positions);

  // Allocate memory for slack variables.
  slack_variables_.reserve(options_.min_num_points_per_view *
                               orientations.size() +
                           view_pairs_.size());

  // Add the constraints to the problem.
  AddCameraToCameraConstraints(orientations, positions);
  if (options_.min_num_points_per_view > 0) {
    AddPointToCameraConstraints(orientations, positions);
  }

  // Set one camera to be at the origin to remove the ambiguity of the origin.
  positions->begin()->second.setZero();
  problem_->SetParameterBlockConstant(positions->begin()->second.data());

  // Set the solver options.
  ceres::Solver::Summary summary;
  solver_options_.logging_type =
      options_.verbose ? ceres::PER_MINIMIZER_ITERATION : ceres::SILENT;
  solver_options_.num_threads = options_.num_threads;
  solver_options_.num_linear_solver_threads = options_.num_threads;
  solver_options_.max_num_iterations = options_.max_num_iterations;

  if (positions->size() > kMinNumCamerasForIterativeSchur) {
    solver_options_.linear_solver_type = ceres::ITERATIVE_SCHUR;
    solver_options_.preconditioner_type = ceres::SCHUR_JACOBI;
  } else {
    solver_options_.linear_solver_type = ceres::SPARSE_SCHUR;
  }

  solver_options_.function_tolerance = 1e-10;
  solver_options_.parameter_tolerance = 1e-10;

  ceres::Solve(solver_options_, problem_.get(), &summary);
  LOG_IF(INFO, options_.verbose) << summary.FullReport();
  if (!summary.IsSolutionUsable()) {
    return false;
  }

  // Reweight the problem and solve again.
  for (int i = 0; i < options_.max_num_reweighted_iterations; i++) {
    ReweightLossFunctions();
    const double cost_change =
        (summary.initial_cost - summary.final_cost) / summary.final_cost;
    if (cost_change < options_.reweighted_convergence_tolerance) {
      VLOG(1) << "IRLS positions solver converged in " << i << " iterations.";
      break;
    }

    // Solve the reweighted problem.
    ceres::Solve(solver_options_, problem_.get(), &summary);
    LOG_IF(INFO, options_.verbose) << summary.FullReport();
    if (!summary.IsSolutionUsable()) {
      return false;
    }
  }

  STLDeleteElements(&slack_variables_);

  return summary.IsSolutionUsable();
}

void NonlinearPositionEstimator::InitializeRandomPositions(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  std::unordered_set<ViewId> constrained_positions;
  constrained_positions.reserve(orientations.size());
  for (const auto& view_pair : view_pairs_) {
    constrained_positions.insert(view_pair.first.first);
    constrained_positions.insert(view_pair.first.second);
  }

  positions->reserve(orientations.size());
  for (const auto& orientation : orientations) {
    if (ContainsKey(constrained_positions, orientation.first)) {
      (*positions)[orientation.first] = 100.0 * Vector3d::Random();
    }
  }
}

void NonlinearPositionEstimator::AddCameraToCameraConstraints(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  for (const auto& view_pair : view_pairs_) {
    const ViewId view_id1 = view_pair.first.first;
    const ViewId view_id2 = view_pair.first.second;
    Vector3d* position1 = FindOrNull(*positions, view_id1);
    Vector3d* position2 = FindOrNull(*positions, view_id2);

    // Do not add this view pair if one or both of the positions do not exist.
    if (position1 == nullptr || position2 == nullptr) {
      continue;
    }

    // Rotate the relative translation so that it is aligned to the global
    // orientation frame.
    const Vector3d translation_direction = GetRotatedTranslation(
        FindOrDie(orientations, view_id1), view_pair.second.position_2);

    double* slack_variable = new double((*position2 - *position1).norm());

    ceres::CostFunction* cost_function =
        PairwiseTranslationError::Create(translation_direction, 1.0);

    // Create loss function wrapper and move it to the back of the loss function
    // container.
    ceres::LossFunctionWrapper* loss_function = new ceres::LossFunctionWrapper(
        new ceres::HuberLoss(options_.robust_loss_width),
        ceres::TAKE_OWNERSHIP);
    loss_functions_.push_back(loss_function);

    problem_->AddResidualBlock(cost_function,
                               loss_functions_.back(),
                               position1->data(),
                               position2->data(),
                               slack_variable);

    // Set lower bound for the slack variables.
    slack_variables_.emplace_back(slack_variable);
    problem_->SetParameterLowerBound(slack_variable, 0, 1.0);
  }

  VLOG(2) << problem_->NumResidualBlocks() << " camera to camera constraints "
                                              "were added to the position "
                                              "estimation problem.";
}

void NonlinearPositionEstimator::AddPointToCameraConstraints(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    std::unordered_map<ViewId, Eigen::Vector3d>* positions) {
  const int num_camera_to_camera_constraints = problem_->NumResidualBlocks();
  std::unordered_set<TrackId> tracks_to_add;
  const int num_point_to_camera_constraints =
      FindTracksForProblem(*positions, &tracks_to_add);
  if (num_point_to_camera_constraints == 0) {
    return;
  }

  const double point_to_camera_weight =
      options_.point_to_camera_weight *
      static_cast<double>(num_camera_to_camera_constraints) /
      static_cast<double>(num_point_to_camera_constraints);

  triangulated_points_.reserve(tracks_to_add.size());
  for (const TrackId track_id : tracks_to_add) {
    triangulated_points_[track_id] = 100.0 * Vector3d::Random();

    AddTrackToProblem(track_id,
                      orientations,
                      point_to_camera_weight,
                      positions);
  }

  VLOG(2) << num_point_to_camera_constraints << " point to camera constriants "
                                                "were added to the position "
                                                "estimation problem.";
}

int NonlinearPositionEstimator::FindTracksForProblem(
    const std::unordered_map<ViewId, Eigen::Vector3d>& positions,
    std::unordered_set<TrackId>* tracks_to_add) {
  CHECK_NOTNULL(tracks_to_add)->clear();

  std::unordered_map<ViewId, int> tracks_per_camera;
  for (const auto& position : positions) {
    tracks_per_camera[position.first] = 0;
  }

  // Add the tracks that see the most views until each camera has the minimum
  // number of tracks.
  for (const auto& position : positions) {
    const View* view = reconstruction_.View(position.first);
    if (view == nullptr ||
        view->NumFeatures() < options_.min_num_points_per_view) {
      continue;
    }

    // Get the tracks in sorted order so that we add the tracks that see the
    // most cameras first.
    const std::vector<TrackId>& sorted_tracks =
        GetTracksSortedByNumViews(reconstruction_, *view, *tracks_to_add);

    for (int i = 0;
         i < sorted_tracks.size() && tracks_per_camera[position.first] <
                                         options_.min_num_points_per_view;
         i++) {
      // Update the number of point to camera constraints for each camera.
      tracks_to_add->insert(sorted_tracks[i]);
      for (const ViewId view_id :
           reconstruction_.Track(sorted_tracks[i])->ViewIds()) {
        if (!ContainsKey(positions, view_id)) {
          continue;
        }
        ++tracks_per_camera[view_id];
      }
    }
  }

  int num_point_to_camera_constraints = 0;
  for (const auto& tracks_in_camera : tracks_per_camera) {
    num_point_to_camera_constraints += tracks_in_camera.second;
  }
  return num_point_to_camera_constraints;
}

std::vector<TrackId> NonlinearPositionEstimator::GetTracksSortedByNumViews(
    const Reconstruction& reconstruction,
    const View& view,
    const std::unordered_set<TrackId>& existing_tracks) {
  std::vector<std::pair<TrackId, int> > views_per_track;
  views_per_track.reserve(view.NumFeatures());
  const auto& track_ids = view.TrackIds();
  for (const auto& track_id : track_ids) {
    const Track* track = reconstruction.Track(track_id);

    if (track == nullptr || ContainsKey(existing_tracks, track_id)) {
      continue;
    }
    views_per_track.emplace_back(track_id, track->NumViews());
  }

  // Return an empty array if no tracks could be found for this view.
  std::vector<TrackId> sorted_tracks(views_per_track.size());
  if (views_per_track.size() == 0) {
    return sorted_tracks;
  }

  // Sort the tracks by the number of views. Only sort the first few tracks
  // since those are the ones that will be added to the problem.
  const int num_tracks_to_sort =
      std::min(static_cast<int>(views_per_track.size()),
               options_.min_num_points_per_view);
  std::partial_sort(views_per_track.begin(),
                    views_per_track.begin() + num_tracks_to_sort,
                    views_per_track.end(), CompareViewsPerTrack);

  for (int i = 0; i < num_tracks_to_sort; i++) {
    sorted_tracks[i] = views_per_track[i].first;
  }
  return sorted_tracks;
}

void NonlinearPositionEstimator::AddTrackToProblem(
    const TrackId track_id,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    const double point_to_camera_weight,
    std::unordered_map<ViewId, Vector3d>* positions) {
  // For each view in the track add the point to camera correspondences.
  for (const ViewId view_id : reconstruction_.Track(track_id)->ViewIds()) {
    if (!ContainsKey(*positions, view_id)) {
      continue;
    }
    Vector3d& camera_position = FindOrDie(*positions, view_id);
    Vector3d& point = FindOrDie(triangulated_points_, track_id);

    // Rotate the feature ray to be in the global orientation frame.
    const Vector3d feature_ray = GetRotatedFeatureRay(
        reconstruction_.View(view_id)->Camera(),
        FindOrDie(orientations, view_id),
        *reconstruction_.View(view_id)->GetFeature(track_id));

    // Create the slack variable.
    double* slack_variable = new double((camera_position - point).norm());

    // Rotate the relative translation so that it is aligned to the global
    // orientation frame.
    ceres::CostFunction* cost_function =
        PairwiseTranslationError::Create(feature_ray, point_to_camera_weight);

    ceres::LossFunctionWrapper* loss_function = new ceres::LossFunctionWrapper(
        new ceres::HuberLoss(options_.robust_loss_width),
        ceres::TAKE_OWNERSHIP);
    loss_functions_.push_back(loss_function);

    // Add the residual block
    problem_->AddResidualBlock(cost_function,
                               loss_functions_.back(),
                               camera_position.data(),
                               point.data(),
                               slack_variable);

    // Set the lower bound for the slack variable.
    slack_variables_.emplace_back(slack_variable);
    problem_->SetParameterLowerBound(slack_variable, 0, 1.0);
  }
}

void NonlinearPositionEstimator::ReweightLossFunctions() {
  // The regularizer prevents any element in the minimization from having an
  // unbounded influence on the problem.
  static const double kRegularizer = 1e-4;

  // Set evaluation options.
  ceres::Problem::EvaluateOptions evaluate_options;
  evaluate_options.apply_loss_function = false;
  evaluate_options.num_threads = options_.num_threads;

  // Evaluate residual block.
  std::vector<double> residuals;
  CHECK(problem_->Evaluate(evaluate_options, NULL, &residuals, NULL, NULL))
      << "Could not evaluate the residuals for IRLS.";

  // Reweight loss function for the residual block.
  const int num_residual_blocks = residuals.size() / 3;
  DCHECK_EQ(num_residual_blocks, loss_functions_.size());
  for (int i = 0; i < num_residual_blocks; i++) {
    const double sq_residual = residuals[3 * i] * residuals[3 * i] +
                            residuals[3 * i + 1] * residuals[3 * i + 1] +
                            residuals[3 * i + 2] * residuals[3 * i + 2];
    const double reweight_scale = 1.0  / std::sqrt(sq_residual + kRegularizer);
    loss_functions_[i]->Reset(
        new ceres::ScaledLoss(NULL, reweight_scale, ceres::TAKE_OWNERSHIP),
        ceres::TAKE_OWNERSHIP);
  }
}

}  // namespace theia
