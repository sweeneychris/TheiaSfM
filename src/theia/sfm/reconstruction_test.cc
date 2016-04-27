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

#include "gtest/gtest.h"

#include "theia/sfm/reconstruction.h"
#include "theia/util/map_util.h"

namespace theia {

const std::vector<std::string> view_names = {"1", "2", "3"};
const std::vector<Feature> features = { Feature(1, 1),
                                        Feature(2, 2),
                                        Feature(3, 3) };

TEST(Reconstruction, ViewIdFromNameValid) {
  Reconstruction reconstruction;
  const ViewId gt_view_id = reconstruction.AddView(view_names[0]);

  const ViewId view_id = reconstruction.ViewIdFromName(view_names[0]);
  EXPECT_EQ(gt_view_id, view_id);
}

TEST(Reconstruction, ViewIdFromNameInvalid) {
  Reconstruction reconstruction;
  EXPECT_EQ(reconstruction.ViewIdFromName(view_names[0]), kInvalidViewId);
}

TEST(Reconstruction, AddView) {
  Reconstruction reconstruction;
  const ViewId view_id = reconstruction.AddView(view_names[0]);
  EXPECT_NE(view_id, kInvalidViewId);
  EXPECT_EQ(reconstruction.NumViews(), 1);
  EXPECT_EQ(reconstruction.NumTracks(), 0);
  EXPECT_EQ(reconstruction.AddView(view_names[0]), kInvalidViewId);
  EXPECT_EQ(reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id), 0);
}

TEST(Reconstruction, AddViewWithCameraIntrinsicsGroup) {
  Reconstruction reconstruction;
  const CameraIntrinsicsGroupId intrinsics_id = 1;
  const ViewId view_id = reconstruction.AddView(view_names[0], intrinsics_id);
  EXPECT_NE(view_id, kInvalidViewId);
  EXPECT_EQ(reconstruction.NumViews(), 1);
  EXPECT_EQ(reconstruction.NumTracks(), 0);
  EXPECT_EQ(reconstruction.NumCameraIntrinsicGroups(), 1);
  EXPECT_EQ(reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id),
            intrinsics_id);
  EXPECT_EQ(reconstruction.AddView(view_names[0]), kInvalidViewId);
}

TEST(Reconstruction, RemoveView) {
  Reconstruction reconstruction;
  const ViewId view_id1 = reconstruction.AddView(view_names[0]);
  const ViewId view_id2 = reconstruction.AddView(view_names[1]);
  EXPECT_EQ(reconstruction.NumViews(), 2);
  EXPECT_EQ(reconstruction.NumCameraIntrinsicGroups(), 2);

  const CameraIntrinsicsGroupId view1_group =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id1);
  const CameraIntrinsicsGroupId view2_group =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id2);

  EXPECT_TRUE(reconstruction.RemoveView(view_id1));
  EXPECT_EQ(reconstruction.NumViews(), 1);
  EXPECT_EQ(reconstruction.ViewIdFromName(view_names[0]), kInvalidViewId);
  EXPECT_EQ(reconstruction.View(view_id1), nullptr);
  EXPECT_EQ(reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id1),
            kInvalidCameraIntrinsicsGroupId);
  EXPECT_EQ(reconstruction.NumCameraIntrinsicGroups(), 1);
  const std::unordered_set<ViewId> view1_camera_intrinsics_group =
      reconstruction.GetViewsInCameraIntrinsicGroup(view1_group);
  EXPECT_FALSE(ContainsKey(view1_camera_intrinsics_group, view_id1));

  EXPECT_TRUE(reconstruction.RemoveView(view_id2));
  EXPECT_EQ(reconstruction.NumViews(), 0);
  EXPECT_EQ(reconstruction.ViewIdFromName(view_names[1]), kInvalidViewId);
  EXPECT_EQ(reconstruction.View(view_id2), nullptr);
  EXPECT_EQ(reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id2),
            kInvalidCameraIntrinsicsGroupId);
  EXPECT_EQ(reconstruction.NumCameraIntrinsicGroups(), 0);
  const std::unordered_set<ViewId> view2_camera_intrinsics_group =
      reconstruction.GetViewsInCameraIntrinsicGroup(view2_group);
  EXPECT_FALSE(ContainsKey(view2_camera_intrinsics_group, view_id2));

  EXPECT_FALSE(reconstruction.RemoveView(kInvalidViewId));
  EXPECT_FALSE(reconstruction.RemoveView(view_id1));
}

TEST(Reconstruction, GetViewValid) {
  Reconstruction reconstruction;
  const ViewId view_id = reconstruction.AddView(view_names[0]);
  EXPECT_NE(view_id, kInvalidViewId);

  const View* const_view = reconstruction.View(view_id);
  EXPECT_NE(const_view, nullptr);

  View* mutable_view = reconstruction.MutableView(view_id);
  EXPECT_NE(mutable_view, nullptr);
}

TEST(Reconstruction, GetViewValidInvalid) {
  Reconstruction reconstruction;
  static const ViewId view_id = 0;
  const View* const_view = reconstruction.View(view_id);
  EXPECT_EQ(const_view, nullptr);

  View* mutable_view = reconstruction.MutableView(view_id);
  EXPECT_EQ(mutable_view, nullptr);
}

TEST(Reconstruction, GetViewsInCameraIntrinsicGroup) {
  Reconstruction reconstruction;
  const ViewId view_id1 = reconstruction.AddView(view_names[0]);
  const CameraIntrinsicsGroupId intrinsics_id1 =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id1);

  // Add a second view with to the same camera intrinsics group.
  const ViewId view_id2 = reconstruction.AddView(view_names[1], intrinsics_id1);
  const CameraIntrinsicsGroupId intrinsics_id2 =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id2);
  EXPECT_EQ(intrinsics_id1, intrinsics_id2);

  // Add a third view that is in it's own camera intrinsics group.
  const ViewId view_id3 = reconstruction.AddView(view_names[2]);
  const CameraIntrinsicsGroupId intrinsics_id3 =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id3);
  EXPECT_NE(intrinsics_id1, intrinsics_id3);
  EXPECT_EQ(reconstruction.NumCameraIntrinsicGroups(), 2);
}

TEST(Reconstruction, CameraIntrinsicsGroupIds) {
  Reconstruction reconstruction;
  const ViewId view_id1 = reconstruction.AddView(view_names[0]);
  const CameraIntrinsicsGroupId intrinsics_id1 =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id1);

  // Add a second view with to the same camera intrinsics group.
  const ViewId view_id2 = reconstruction.AddView(view_names[1], intrinsics_id1);
  const CameraIntrinsicsGroupId intrinsics_id2 =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id2);
  EXPECT_EQ(intrinsics_id1, intrinsics_id2);

  // Add a third view that is in it's own camera intrinsics group.
  const ViewId view_id3 = reconstruction.AddView(view_names[2]);
  const CameraIntrinsicsGroupId intrinsics_id3 =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id3);
  EXPECT_NE(intrinsics_id1, intrinsics_id3);
  EXPECT_EQ(reconstruction.NumCameraIntrinsicGroups(), 2);

  // Ensure that the group ids are correct.
  const std::unordered_set<CameraIntrinsicsGroupId> group_ids =
      reconstruction.CameraIntrinsicsGroupIds();
  EXPECT_EQ(group_ids.size(), 2);
  EXPECT_TRUE(ContainsKey(group_ids, intrinsics_id1));
  EXPECT_TRUE(ContainsKey(group_ids, intrinsics_id3));
}

TEST(Reconstruction, AddTrackValid) {
  Reconstruction reconstruction;

  const std::vector<std::pair<ViewId, Feature> > track = {
    { 0, features[0] }, { 1, features[1] }
  };
  EXPECT_NE(reconstruction.AddView(view_names[0]), kInvalidViewId);
  EXPECT_NE(reconstruction.AddView(view_names[1]), kInvalidViewId);

  const TrackId track_id = reconstruction.AddTrack(track);
  EXPECT_NE(track_id, kInvalidTrackId);
  EXPECT_TRUE(reconstruction.Track(track_id) != nullptr);
  EXPECT_EQ(reconstruction.NumTracks(), 1);
}

TEST(Reconstruction, AddTrackInvalid) {
  Reconstruction reconstruction;

  // Should fail with less than two views.
  const std::vector<std::pair<ViewId, Feature> > small_track = {
    { 0, features[0] }
  };
  EXPECT_NE(reconstruction.AddView(view_names[0]), kInvalidViewId);
  EXPECT_EQ(reconstruction.AddTrack(small_track), kInvalidTrackId);
  EXPECT_EQ(reconstruction.NumTracks(), 0);
}

TEST(Reconstruction, RemoveTrackValid) {
  Reconstruction reconstruction;

  const std::vector<std::pair<ViewId, Feature> > track = {
    { 0, features[0] }, { 1, features[1] }
  };

  // Should be able to successfully remove the track.
  EXPECT_NE(reconstruction.AddView(view_names[0]), kInvalidViewId);
  EXPECT_NE(reconstruction.AddView(view_names[1]), kInvalidViewId);
  const TrackId track_id = reconstruction.AddTrack(track);
  EXPECT_TRUE(reconstruction.RemoveTrack(track_id));
}

TEST(Reconstruction, RemoveTrackInvalid) {
  Reconstruction reconstruction;

  // Should return false when trying to remove a track not in the
  // reconstruction.
  EXPECT_FALSE(reconstruction.RemoveTrack(kInvalidTrackId));
}

TEST(Reconstruction, GetTrackValid) {
  Reconstruction reconstruction;
  const std::vector<std::pair<ViewId, Feature> > track = {
    { 0, features[0] }, { 1, features[1] }
  };
  EXPECT_NE(reconstruction.AddView(view_names[0]), kInvalidViewId);
  EXPECT_NE(reconstruction.AddView(view_names[1]), kInvalidViewId);
  const TrackId track_id = reconstruction.AddTrack(track);
  EXPECT_NE(track_id, kInvalidTrackId);

  const Track* const_track = reconstruction.Track(track_id);
  EXPECT_NE(const_track, nullptr);

  Track* mutable_track = reconstruction.MutableTrack(track_id);
  EXPECT_NE(mutable_track, nullptr);
}

TEST(Reconstruction, GetTrackInvalid) {
  Reconstruction reconstruction;
  const std::vector<std::pair<ViewId, Feature> > track = {};
  const TrackId track_id = reconstruction.AddTrack(track);
  EXPECT_EQ(track_id, kInvalidTrackId);

  const Track* const_track = reconstruction.Track(track_id);
  EXPECT_EQ(const_track, nullptr);

  Track* mutable_track = reconstruction.MutableTrack(track_id);
  EXPECT_EQ(mutable_track, nullptr);
}

}  // namespace theia
