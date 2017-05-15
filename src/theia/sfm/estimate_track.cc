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

#include "theia/sfm/estimate_track.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
#include <vector>

#include "theia/math/util.h"
#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/estimators/estimate_triangulation.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/track.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/sfm/types.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"

namespace theia {

namespace {

void GetObservationsFromTrackViews(
    const TrackId track_id,
    const Reconstruction& reconstruction,
    std::vector<ViewId>* view_ids,
    std::vector<Eigen::Vector2d>* features,
    std::vector<Eigen::Vector3d>* origins,
    std::vector<Eigen::Vector3d>* ray_directions) {
  const Track track = *reconstruction.Track(track_id);
  for (const ViewId view_id : track.ViewIds()) {
    const View* view = reconstruction.View(view_id);

    // Skip this view if it does not exist or has not been estimated yet.
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }

    // If the feature is not in the view then we have an ill-formed
    // reconstruction.
    const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
    const Eigen::Vector3d image_ray =
        view->Camera().PixelToUnitDepthRay(*feature).normalized();

    features->emplace_back(*feature);
    view_ids->emplace_back(view_id);
    origins->emplace_back(view->Camera().GetPosition());
    ray_directions->emplace_back(image_ray);
  }
}

// Returns false if the reprojection error of the triangulated point is greater
// than the max allowable reprojection error (for any observation) and true
// otherwise.
bool AcceptableReprojectionError(
    const Reconstruction& reconstruction,
    const TrackId& track_id,
    const std::vector<ViewId>& view_ids,
    const std::vector<Eigen::Vector2d>& features,
    const double sq_max_reprojection_error_pixels) {
  const Track& track = *reconstruction.Track(track_id);
  int num_projections = 0;
  double mean_sq_reprojection_error = 0;
  for (int i = 0; i < view_ids.size(); i++) {
    // We do not need to check if the view is present or estimated, since this
    // was done prior to the input.
    const View* view = reconstruction.View(view_ids[i]);
    const Camera& camera = view->Camera();
    Eigen::Vector2d reprojection;
    if (camera.ProjectPoint(track.Point(), &reprojection) < 0) {
      return false;
    }

    mean_sq_reprojection_error += (features[i] - reprojection).squaredNorm();
    ++num_projections;
  }

  return (mean_sq_reprojection_error / static_cast<double>(num_projections)) <
         sq_max_reprojection_error_pixels;
}

}  // namespace

// Estimate only the tracks supplied by the user.
TrackEstimator::Summary TrackEstimator::EstimateAllTracks() {
  // Estimate all tracks that are seen by estimated views.
  const auto& view_ids = reconstruction_->ViewIds();
  std::unordered_set<TrackId> tracks;
  for (const ViewId view_id : view_ids) {
    const View* view = reconstruction_->View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }
    const auto& tracks_in_view = view->TrackIds();
    tracks.insert(tracks_in_view.begin(), tracks_in_view.end());
  }

  return EstimateTracks(tracks);
}

TrackEstimator::Summary TrackEstimator::EstimateTracks(
    const std::unordered_set<TrackId>& track_ids) {
  tracks_to_estimate_.clear();
  summary_ = TrackEstimator::Summary();

  // Get all unestimated track ids.
  tracks_to_estimate_.reserve(track_ids.size());
  for (const TrackId track_id : track_ids) {
    Track* track = reconstruction_->MutableTrack(track_id);
    if (!track->IsEstimated()) {
      tracks_to_estimate_.emplace_back(track_id);
    }
  }
  summary_.input_num_estimated_tracks =
      track_ids.size() - tracks_to_estimate_.size();
  summary_.num_triangulation_attempts = tracks_to_estimate_.size();

  // Exit early if there are no tracks to estimate.
  if (tracks_to_estimate_.size() == 0) {
    return summary_;
  }

  // Estimate the tracks in parallel. Instead of 1 threadpool worker per track,
  // we let each worker estimate a fixed number of tracks at a time (e.g. 20
  // tracks). Since estimating the tracks is so fast, this strategy is better
  // helps speed up multithreaded estimation by reducing the overhead of
  // starting/stopping threads.
  const int num_threads = std::min(
      options_.num_threads, static_cast<int>(tracks_to_estimate_.size()));
  const int interval_step =
      std::min(options_.multithreaded_step_size,
               static_cast<int>(tracks_to_estimate_.size()) / num_threads);

  std::unique_ptr<ThreadPool> pool(new ThreadPool(num_threads));
  for (int i = 0; i < tracks_to_estimate_.size(); i += interval_step) {
    const int end_interval = std::min(
        static_cast<int>(tracks_to_estimate_.size()), i + interval_step);
    pool->Add(&TrackEstimator::EstimateTrackSet, this, i, end_interval);
  }

  // Wait for all tracks to be estimated.
  pool.reset(nullptr);

  LOG(INFO) << summary_.estimated_tracks.size() << " tracks were estimated of "
            << summary_.num_triangulation_attempts << " possible tracks.";
  return summary_;
}

void TrackEstimator::EstimateTrackSet(const int start, const int end) {
  std::unordered_set<TrackId> estimated_tracks;
  for (int i = start; i < end; i++) {
    if (EstimateTrack(tracks_to_estimate_[i])) {
      estimated_tracks.emplace(tracks_to_estimate_[i]);
    }
  }

  // Add the estimated tracks to the output summary.
  std::lock_guard<std::mutex> guard(summary_mutex_);
  summary_.estimated_tracks.insert(estimated_tracks.begin(),
                                   estimated_tracks.end());
}

bool TrackEstimator::EstimateTrack(const TrackId track_id) {
  static const int kMinNumObservationsForTriangulation = 2;

  Track* track = reconstruction_->MutableTrack(track_id);
  CHECK(!track->IsEstimated()) << "Track " << track_id
                               << " is already estimated.";

  // Gather projection matrices and features.
  std::vector<ViewId> view_ids;
  std::vector<Eigen::Vector2d> features;
  std::vector<Eigen::Vector3d> origins, ray_directions;
  GetObservationsFromTrackViews(track_id,
                                *reconstruction_,
                                &view_ids,
                                &features,
                                &origins,
                                &ray_directions);

  // Check the angle between views.
  if (view_ids.size() < kMinNumObservationsForTriangulation ||
      !SufficientTriangulationAngle(ray_directions,
                                    options_.min_triangulation_angle_degrees)) {
    return false;
  }

  // Triangulate the track.
  if (!TriangulateMidpoint(origins, ray_directions, track->MutablePoint())) {
    return false;
  }

  // Bundle adjust the track.
  if (options_.bundle_adjustment) {
    track->SetEstimated(true);
    const BundleAdjustmentSummary summary =
        BundleAdjustTrack(options_.ba_options, track_id, reconstruction_);
    track->SetEstimated(false);
    if (!summary.success) {
      return false;
    }
  }

  // Ensure the reprojection errors are acceptable.
  const double sq_max_reprojection_error_pixels =
      options_.max_acceptable_reprojection_error_pixels *
      options_.max_acceptable_reprojection_error_pixels;

  if (!AcceptableReprojectionError(*reconstruction_,
                                   track_id,
                                   view_ids,
                                   features,
                                   sq_max_reprojection_error_pixels)) {
    return false;
  }

  track->SetEstimated(true);
  return true;
}

}  // namespace theia
