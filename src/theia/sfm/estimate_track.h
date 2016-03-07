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

#ifndef THEIA_SFM_ESTIMATE_TRACK_H_
#define THEIA_SFM_ESTIMATE_TRACK_H_

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"

#include <unordered_set>

namespace theia {

// Estimates the 3D point of a track by using all estimated views to compute a
// (potentially nonminimal) triangulation of track. The the angle between all
// views and the triangulated point must be greater than the minimum
// triangulation error. The track estimation is successful if all views have a
// reprojection error less than the specified max reprojection error.
// Estimates all unestimated tracks in the reconstruction.
class TrackEstimator {
 public:
  struct Options {
    // Number of threads for multithreading.
    int num_threads = 1;

    // Maximum reprojection error for successful triangulation.
    double max_acceptable_reprojection_error_pixels = 5.0;

    // Minimum triangulation angle between two views required for
    // triangulation. For N-View triangulation we require that at least one pair
    // of views has this an angle this large.
    double min_triangulation_angle_degrees = 3.0;

    // Perform bundle adjustment on the track as soon as a position is
    // estimated.
    bool bundle_adjustment = true;
    BundleAdjustmentOptions ba_options;

    // For thread-level parallelism, it is better to estimate a small fixed
    // number of tracks per thread worker instead of 1 track per worker. This
    // number controls how many points are estimated per worker.
    int multithreaded_step_size = 100;
  };

  struct Summary {
    // Number of estimated tracks that were input.
    int input_num_estimated_tracks = 0;

    // Number of estimated tracks after trianguation.
    int final_num_estimated_tracks = 0;

    // Number of triangulation attempts made.
    int num_triangulation_attempts = 0;

    // TrackId of the newly estimated tracks. This set does not include tracks
    // that were input as estimated.
    std::unordered_set<TrackId> estimated_tracks;

    double ba_setup_time_in_seconds = 0;
    double ba_total_time_in_seconds = 0;
  };

  TrackEstimator(const Options& options, Reconstruction* reconstruction)
      : options_(options), reconstruction_(reconstruction) {}

  // Attempts to estimate all unestimated tracks.
  Summary EstimateAllTracks();

  // Estimate only the tracks supplied by the user.
  Summary EstimateTracks(const std::unordered_set<TrackId>& track_ids);

 private:
  void EstimateTrackSet(const int start, const int stop);
  bool EstimateTrack(const TrackId track_id);

  const Options options_;
  Reconstruction* reconstruction_;
  std::vector<TrackId> tracks_to_estimate_;
  std::unordered_set<TrackId> successfully_estimated_tracks_;
};

}  // namespace theia

#endif  // THEIA_SFM_ESTIMATE_TRACK_H_
