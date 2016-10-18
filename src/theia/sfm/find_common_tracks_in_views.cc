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

#include "theia/sfm/find_common_tracks_in_views.h"

#include <glog/logging.h>
#include <iterator>
#include <algorithm>
#include <vector>

#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {
namespace {
std::vector<TrackId> GetTrackIdsForView(const Reconstruction& reconstruction,
                                        const ViewId view_id) {
  const View* view = CHECK_NOTNULL(reconstruction.View(view_id));
  return view->TrackIds();
}

}  // namespace

// Finds the tracks that are common to all views and returns them. An empty
// vector is returned if no tracks are observed in all of the views.
std::vector<TrackId> FindCommonTracksInViews(
    const Reconstruction& reconstruction, const std::vector<ViewId>& views) {
  CHECK_GT(views.size(), 1)
      << "Finding common tracks between views requires at least 2 views.";

  std::vector<TrackId> common_track_ids =
      GetTrackIdsForView(reconstruction, views[0]);
  std::sort(common_track_ids.begin(), common_track_ids.end());

  std::vector<TrackId> buffer;
  for (int i = 1; i < views.size(); i++) {
    if (common_track_ids.size() == 0) {
      break;
    }

    buffer.clear();

    std::vector<TrackId> sorted_track_ids =
        GetTrackIdsForView(reconstruction, views[i]);
    std::sort(sorted_track_ids.begin(), sorted_track_ids.end());

    std::set_intersection(common_track_ids.begin(), common_track_ids.end(),
                          sorted_track_ids.begin(), sorted_track_ids.end(),
                          std::back_inserter(buffer));
    std::swap(common_track_ids, buffer);
  }

  return common_track_ids;
}

}  // namespace theia
