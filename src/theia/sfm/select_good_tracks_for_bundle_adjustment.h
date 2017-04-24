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

#ifndef THEIA_SFM_CHOOSE_GOOD_TRACKS_FOR_BUNDLE_ADJUSTMENT_H_
#define THEIA_SFM_SELECT_GOOD_TRACKS_FOR_BUNDLE_ADJUSTMENT_H_

#include <unordered_set>

#include "theia/sfm/types.h"

namespace theia {
class Reconstruction;

// The efficiency of large scale bundle adjustment can be dramatically increased
// by choosing only a subset of 3d points to optimize, as the 3d points tend to
// have increasing scene redundancy. If the points are chosen in a way that
// properly constrains the nonlinear optimization, similar results in quality
// may be observed compared to optimized all tracks.
//
// The 3d points are chosen such that they fit the following criteria:
//    a) High confidence (i.e. low reprojection error).
//    b) Long tracks are preferred.
//    c) The tracks used for optimization provide a good spatial coverage in
//       each image.
//    d) Each view observes at least K optimized tracks.
//
// Tracks in each image are first hashed into spatial bins with an image grid
// where each image grid cell is the provided width. Within each grid cell, the
// tracks are ordered based on their track length, then by mean reprojection
// error. The track length is truncated to be no longer than
// long_track_length_threshold so that among long tracks, the ones with low
// reprojection error are chosen for bundle adjustment.
//
// We recommend the grid cell size is set to 100 pixels, the long track length
// threshold is set to 10, and the min num optimized tracks per view is set to
// 100.
bool SelectGoodTracksForBundleAdjustment(
    const Reconstruction& reconstruction,
    const int long_track_length_threshold,
    const int image_grid_cell_size_pixels,
    const int min_num_optimized_tracks_per_view,
    std::unordered_set<TrackId>* tracks_to_optimize);

}  // namespace theia

#endif  // THEIA_SFM_SELECT_GOOD_TRACKS_FOR_BUNDLE_ADJUSTMENT_H_
