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

#include <algorithm>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "theia/sfm/find_common_tracks_in_views.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"

namespace theia {

TEST(FindCommonTracksInViews, NoCommonTracks) {
  const Feature feature;
  const std::vector<std::pair<ViewId, Feature> > track1 = {{0, feature},
                                                           {1, feature}};
  const std::vector<std::pair<ViewId, Feature> > track2 = {{1, feature},
                                                           {2, feature}};
  const std::vector<std::pair<ViewId, Feature> > track3 = {{2, feature},
                                                           {3, feature}};

  Reconstruction reconstruction;
  reconstruction.AddView("0");
  reconstruction.AddView("1");
  reconstruction.AddView("2");
  reconstruction.AddView("3");
  reconstruction.AddTrack(track1);
  reconstruction.AddTrack(track2);
  reconstruction.AddTrack(track3);

  const std::vector<ViewId> views_to_find_common_tracks = { 0, 1, 2, 3};
  const std::vector<TrackId> common_tracks = FindCommonTracksInViews(
      reconstruction, views_to_find_common_tracks);
  EXPECT_EQ(common_tracks.size(), 0);
}


TEST(FindCommonTracksInViews, CommonTracks) {
  const Feature feature;
  const std::vector<std::pair<ViewId, Feature> > track1 = {{1, feature},
                                                           {2, feature}};
  const std::vector<std::pair<ViewId, Feature> > track2 = {
      {0, feature}, {1, feature}, {2, feature}, {3, feature}};
  const std::vector<std::pair<ViewId, Feature> > track3 = {{2, feature},
                                                           {3, feature}};

  Reconstruction reconstruction;
  reconstruction.AddView("0");
  reconstruction.AddView("1");
  reconstruction.AddView("2");
  reconstruction.AddView("3");
  reconstruction.AddTrack(track1);
  reconstruction.AddTrack(track2);
  reconstruction.AddTrack(track3);

  const std::vector<ViewId> views_to_find_common_tracks = { 0, 1, 2, 3};
  const std::vector<TrackId> common_tracks = FindCommonTracksInViews(
      reconstruction, views_to_find_common_tracks);
  EXPECT_EQ(common_tracks.size(), 1);
}

}  // namespace theia
