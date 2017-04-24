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
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include "theia/sfm/set_outlier_tracks_to_unestimated.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <unordered_set>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"
#include "theia/sfm/triangulation/triangulation.h"

namespace theia {

int SetOutlierTracksToUnestimated(const double max_inlier_reprojection_error,
                                  const double min_triangulation_angle_degrees,
                                  Reconstruction* reconstruction) {
  const auto& track_ids = reconstruction->TrackIds();
  const std::unordered_set<TrackId> all_tracks(track_ids.begin(),
                                               track_ids.end());
  return SetOutlierTracksToUnestimated(all_tracks,
                                       max_inlier_reprojection_error,
                                       min_triangulation_angle_degrees,
                                       reconstruction);
}

int SetOutlierTracksToUnestimated(const std::unordered_set<TrackId>& track_ids,
                                  const double max_inlier_reprojection_error,
                                  const double min_triangulation_angle_degrees,
                                  Reconstruction* reconstruction) {
  const double max_sq_reprojection_error =
      max_inlier_reprojection_error * max_inlier_reprojection_error;

  int num_estimated_tracks = 0;
  int num_bad_reprojections = 0;
  int num_insufficient_viewing_angles = 0;

  for (const TrackId track_id : track_ids) {
    Track* track = reconstruction->MutableTrack(track_id);
    if (!track->IsEstimated()) {
      continue;
    }
    ++num_estimated_tracks;

    std::vector<Eigen::Vector3d> ray_directions;
    const auto& view_ids = track->ViewIds();
    int num_projections = 0;
    double mean_sq_reprojection_error = 0;
    for (const ViewId view_id : view_ids) {
      const View* view = CHECK_NOTNULL(reconstruction->View(view_id));
      if (!view->IsEstimated()) {
        continue;
      }

      const Camera& camera = view->Camera();
      const Eigen::Vector3d ray_direction =
          track->Point().hnormalized() - camera.GetPosition();
      ray_directions.push_back(ray_direction.normalized());

      // Check the reprojection error.
      const Feature* feature = view->GetFeature(track_id);
      // Reproject the observations.
      Eigen::Vector2d projection;
      const double depth = camera.ProjectPoint(track->Point(), &projection);
      // Remove the feature if the reprojection is behind the camera.
      if (depth < 0) {
        ++num_bad_reprojections;
        track->SetEstimated(false);
        break;
      }
      mean_sq_reprojection_error += (projection - *feature).squaredNorm();
      ++num_projections;
    }

    mean_sq_reprojection_error /= static_cast<double>(num_projections);
    if (track->IsEstimated() &&
        mean_sq_reprojection_error > max_sq_reprojection_error) {
      ++num_bad_reprojections;
      track->SetEstimated(false);
    }

    // The track will remain estimated if the reprojection errors were all
    // good. We then test that the track is properly constrained by having at
    // least two cameras view it with a sufficient viewing angle.
    if (track->IsEstimated() &&
        !SufficientTriangulationAngle(ray_directions,
                                      min_triangulation_angle_degrees)) {
      ++num_insufficient_viewing_angles;
      track->SetEstimated(false);
    }
  }

  LOG_IF(INFO, num_bad_reprojections > 0 || num_insufficient_viewing_angles > 0)
      << num_bad_reprojections
      << " points were removed because of bad reprojection errors. "
      << num_insufficient_viewing_angles
      << " points were removed because they had insufficient viewing angles "
         "and were poorly constrained.";

  return num_bad_reprojections + num_insufficient_viewing_angles;
}

}  // namespace theia
