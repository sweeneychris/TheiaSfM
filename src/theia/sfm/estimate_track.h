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

#include "theia/sfm/types.h"

namespace theia {

class Reconstruction;

struct EstimateTrackOptions {
  // Maximum reprojection error for successful triangulation.
  double max_acceptable_reprojection_error_pixels = 5.0;

  // Minimum triangulation angle between two views required for
  // triangulation. For N-View triangulation we require that at least one pair
  // of views has this an angle this large.
  double min_triangulation_angle_degrees = 3.0;

  // Perform bundle adjustment on the track as soon as a position is estimated.
  bool bundle_adjustment = true;
};

// Estimates the 3D point of a track by using all estimated views to compute a
// (potentially nonminimal) triangulation of track. The the angle between all
// views and the triangulated point must be greater than the minimum
// triangulation error. The track estimation is successful if all views have a
// reprojection error less than the specified max reprojection error.
bool EstimateTrack(const EstimateTrackOptions& options,
                   const TrackId track_id,
                   Reconstruction* reconstruction);

// Estimates all unestimated tracks in the reconstruction.
void EstimateAllTracks(const EstimateTrackOptions& options,
                       const int num_threads,
                       Reconstruction* reconstruction);


}  // namespace theia

#endif  // THEIA_SFM_ESTIMATE_TRACK_H_
