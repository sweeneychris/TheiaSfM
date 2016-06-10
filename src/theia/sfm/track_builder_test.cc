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

#include <glog/logging.h>

#include <string>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"

#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/track_builder.h"
#include "theia/sfm/types.h"

namespace theia {
static const int kMinTrackLength = 2;

// Ensure that each track has been added to every view.
void VerifyTracks(const Reconstruction& reconstruction) {
  const std::vector<TrackId> track_ids = reconstruction.TrackIds();
  for (const TrackId track_id : track_ids) {
    const Track* track = CHECK_NOTNULL(reconstruction.Track(track_id));
    for (const ViewId view_id : track->ViewIds()) {
      const View* view = CHECK_NOTNULL(reconstruction.View(view_id));
      EXPECT_NE(view->GetFeature(track_id), nullptr);
    }
  }
}

// Perfect tracks.
TEST(TrackBuilder, ConsistentTracks) {
  static const int kMaxTrackLength = 10;
  static const int kNumCorrespondences = 4;

  const ViewId view_ids[kNumCorrespondences][2] = {
    { 0, 1 }, { 0, 1 }, { 1, 2 }, { 1, 2 }
  };

  TrackBuilder track_builder(kMinTrackLength, kMaxTrackLength);
  for (int i = 0; i < kNumCorrespondences; i++) {
    track_builder.AddFeatureCorrespondence(view_ids[i][0],
                                           Feature(i, i),
                                           view_ids[i][1],
                                           Feature(i, i));
  }

  Reconstruction reconstruction;
  reconstruction.AddView("0");
  reconstruction.AddView("1");
  reconstruction.AddView("2");
  track_builder.BuildTracks(&reconstruction);
  VerifyTracks(reconstruction);
  EXPECT_EQ(reconstruction.NumTracks(), kNumCorrespondences);
}

// Singleton tracks.
TEST(TrackBuilder, SingletonTracks) {
  // Having a small max track length will force a singleton track.
  static const int kMaxTrackLength = 2;
  static const int kNumCorrespondences = 2;

  const ViewId view_ids[kNumCorrespondences][2] = {
    { 0, 1 }, { 1, 2 } };

  TrackBuilder track_builder(kMinTrackLength, kMaxTrackLength);
  for (int i = 0; i < kNumCorrespondences; i++) {
    track_builder.AddFeatureCorrespondence(
        view_ids[i][0], Feature(0, 0),
        view_ids[i][1], Feature(0, 0));
  }

  Reconstruction reconstruction;
  reconstruction.AddView("0");
  reconstruction.AddView("1");
  reconstruction.AddView("2");

  track_builder.BuildTracks(&reconstruction);
  VerifyTracks(reconstruction);
  EXPECT_EQ(reconstruction.NumTracks(), 1);
}

// Inconsistent tracks.
TEST(TrackBuilder, InconsistentTracks) {
  static const int kMaxTrackLength = 10;
  static const int kNumCorrespondences = 4;

  const ViewId view_ids[kNumCorrespondences][2] = {
    { 0, 1 }, { 0, 1 }, { 1, 2 }, { 1, 2 }
  };

  TrackBuilder track_builder(kMinTrackLength, kMaxTrackLength);
  for (int i = 0; i < kNumCorrespondences; i++) {
    track_builder.AddFeatureCorrespondence(
        view_ids[i][0], Feature(0, 0),
        view_ids[i][1], Feature(i + 1, i + 1));
  }

  Reconstruction reconstruction;
  reconstruction.AddView("0");
  reconstruction.AddView("1");
  reconstruction.AddView("2");

  track_builder.BuildTracks(&reconstruction);
  VerifyTracks(reconstruction);
  EXPECT_EQ(reconstruction.NumTracks(), 2);
}

// Tracks limited by size.
TEST(TrackBuilder, MaxTrackLength) {
  static const int kMaxTrackLength = 2;
  static const int kNumViews = 6;

  const ViewId view_ids[kNumViews] = { 0, 1, 2, 3, 4, 5 };

  TrackBuilder track_builder(kMinTrackLength, kMaxTrackLength);
  for (int i = 0; i < kNumViews - 1; i++) {
    track_builder.AddFeatureCorrespondence(
        view_ids[i], Feature(0, 0),
        view_ids[i + 1], Feature(0, 0));
  }

  Reconstruction reconstruction;
  reconstruction.AddView("0");
  reconstruction.AddView("1");
  reconstruction.AddView("2");
  reconstruction.AddView("3");
  reconstruction.AddView("4");
  reconstruction.AddView("5");

  track_builder.BuildTracks(&reconstruction);
  VerifyTracks(reconstruction);
  EXPECT_EQ(reconstruction.NumTracks(), 3);
}

TEST(TrackBuilder, MinTrackLength) {
  static const int kMaxTrackLength = 10;
  static const int min_track_length = 3;
  static const int kNumViews = 6;

  const ViewId view_ids[kNumViews] = { 0, 1, 2, 3, 4, 5 };

  TrackBuilder track_builder(min_track_length, kMaxTrackLength);

  // Add one track that is larger than the min track length.
  for (int i = 0; i < kNumViews - 1; i++) {
    track_builder.AddFeatureCorrespondence(
        view_ids[i], Feature(0, 0),
        view_ids[i + 1], Feature(0, 0));
  }

  // Add another track that is smaller than the min track length.
  track_builder.AddFeatureCorrespondence(view_ids[0], Feature(1, 1),
                                         view_ids[1], Feature(1, 1));

  Reconstruction reconstruction;
  reconstruction.AddView("0");
  reconstruction.AddView("1");
  reconstruction.AddView("2");
  reconstruction.AddView("3");
  reconstruction.AddView("4");
  reconstruction.AddView("5");

  track_builder.BuildTracks(&reconstruction);
  VerifyTracks(reconstruction);
  EXPECT_EQ(reconstruction.NumTracks(), 1);
}

}  // namespace theia
