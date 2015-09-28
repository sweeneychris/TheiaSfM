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
#include <vector>

#include "theia/math/util.h"
#include "theia/sfm/bundle_adjustment/bundle_adjust_track.h"
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
    const double sq_max_reprojection_error_pixels) {
  const Track& track = *reconstruction.Track(track_id);
  for (const ViewId view_id : track.ViewIds()) {
    const View* view = reconstruction.View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }
    const Camera& camera = view->Camera();
    const Feature* feature = view->GetFeature(track_id);
    Eigen::Vector2d reprojection;
    if (camera.ProjectPoint(track.Point(), &reprojection) < 0) {
      return false;
    }
    const double sq_reprojection_error =
        (*feature - reprojection).squaredNorm();
    if (sq_reprojection_error > sq_max_reprojection_error_pixels) {
      return false;
    }
  }
  return true;
}

int NumEstimatedViewsObservingTrack(const Reconstruction& reconstruction,
                                    const Track& track) {
  int num_estimated_views = 0;
  for (const ViewId view_id : track.ViewIds()) {
    const View* view = reconstruction.View(view_id);
    if (view != nullptr && view->IsEstimated()) {
      ++num_estimated_views;
    }
  }
  return num_estimated_views;
}

}  // namespace

bool EstimateTrack(const EstimateTrackOptions& options,
                   const TrackId track_id,
                   Reconstruction* reconstruction) {
  Track* track = reconstruction->MutableTrack(track_id);
  CHECK(!track->IsEstimated()) << "Track " << track_id
                               << " is already estimated.";

  // Gather projection matrices and features.
  std::vector<ViewId> view_ids;
  std::vector<Eigen::Vector3d> origins, ray_directions;
  GetObservationsFromTrackViews(track_id,
                                *reconstruction,
                                &view_ids,
                                &origins,
                                &ray_directions);

  // Triangulate with 2 or n views.
  if (origins.size() < 2) {
    VLOG(3) << "Triangulation of track " << track_id
            << " requires at least 2 observations from estimated views.";
    return false;
  }

  // Check the angle between views.
  if (!SufficientTriangulationAngle(ray_directions,
                                    options.min_triangulation_angle_degrees)) {
    VLOG(3)
        << "Track " << track_id
        << " has an insufficient triangulation angle and cannot be estimated.";
    return false;
  }

  // Triangulate the track.
  if (!TriangulateMidpoint(origins, ray_directions, track->MutablePoint())) {
    return false;
  }

  // Bundle adjust the track. The 2-view triangulation method is optimal so we
  // do not need to run BA for that case.
  if (options.bundle_adjustment) {
    track->SetEstimated(true);
    const bool ba_success = BundleAdjustTrack(track_id, reconstruction);
    track->SetEstimated(false);
    if (!ba_success) {
      return false;
    }
  }

  // Ensure the reprojection errors are acceptable.
  const double sq_max_reprojection_error_pixels =
      options.max_acceptable_reprojection_error_pixels *
      options.max_acceptable_reprojection_error_pixels;

  if (!AcceptableReprojectionError(*reconstruction,
                                   track_id,
                                   sq_max_reprojection_error_pixels)) {
    return false;
  }

  track->SetEstimated(true);
  return true;
}

void EstimateAllTracks(const EstimateTrackOptions& options,
                       const int num_threads,
                       Reconstruction* reconstruction) {
  int num_triangulation_attempts = 0;
  std::unique_ptr<ThreadPool> pool(new ThreadPool(num_threads));

  const auto& track_ids = reconstruction->TrackIds();
  int num_estimated_tracks_before = 0;
  for (const TrackId track_id : track_ids) {
    Track* track = reconstruction->MutableTrack(track_id);
    if (track->IsEstimated()) {
      ++num_estimated_tracks_before;
      continue;
    }

    const int num_views_observing_track =
        NumEstimatedViewsObservingTrack(*reconstruction, *track);
    // Skip tracks that do not have enough observations.
    if (num_views_observing_track < 2) {
      continue;
    }

    ++num_triangulation_attempts;
    pool->Add(EstimateTrack, options, track_id, reconstruction);
  }
  pool.reset(nullptr);

  const int num_estimated_tracks =
      NumEstimatedTracks(*reconstruction) - num_estimated_tracks_before;
  LOG(INFO) << num_estimated_tracks << " tracks were estimated of "
            << num_triangulation_attempts << " possible tracks.";
}

}  // namespace theia
