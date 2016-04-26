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

#ifndef THEIA_SFM_RECONSTRUCTION_H_
#define THEIA_SFM_RECONSTRUCTION_H_

#include <cereal/access.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {

// A Reconstruction is the main object for structure from motion. It contains
// all the 3D and camera information as well as visibility constraints. That is,
// it contains Tracks and Views and the constraints between them.
//
// The API of this class is based on LibMV's Reconstruction:
// https://github.com/libmv/libmv/blob/master/src/libmv/reconstruction/reconstruction.h
//
// The main difference is that LibMV keeps the relationship between tracks and
// views in one large global map, whereas we maintain the view to track
// relationship within each view or track object.
class Reconstruction {
 public:
  Reconstruction();
  ~Reconstruction();

  // Returns the unique ViewId of the view name, or kInvalidViewId if the view
  // does not
  // exist.
  ViewId ViewIdFromName(const std::string& view_name) const;

  // Creates a new view and returns the view id.
  ViewId AddView(const std::string& view_name);
  // Creates a new view and assigns it to the specified camera intrinsics group.
  ViewId AddView(const std::string& view_name,
                 const CameraIntrinsicsGroupId group_id);

  // Removes the view from the reconstruction and removes all references to the
  // view in the tracks. Any tracks that have zero views after this view is
  // removed are also removed.
  bool RemoveView(const ViewId view_id);
  int NumViews() const;

  // Returns the View or a nullptr if the track does not exist.
  const class View* View(const ViewId view_id) const;
  class View* MutableView(const ViewId view_id);

  // Return all ViewIds in the reconstruction.
  std::vector<ViewId> ViewIds() const;

  // Get the camera intrinsics group id for the view id.
  CameraIntrinsicsGroupId CameraIntrinsicsGroupIdFromViewId(
      const ViewId view_id) const;

  // Return all view ids with the given camera intrinsics group id. If an
  // invalid or non-existant group is chosen then an empty set will be returned.
  std::unordered_set<ViewId> GetViewsInCameraIntrinsicGroup(
      const CameraIntrinsicsGroupId group_id) const;
  int NumCameraIntrinsicGroups() const;

  // Returns all group ids.
  std::unordered_set<CameraIntrinsicsGroupId> CameraIntrinsicsGroupIds() const;

  // Add a new track to the reconstruction. If successful, the new track id is
  // returned. Failure results when multiple features from the same image are
  // present, and kInvalidTrackId is returned.
  TrackId AddTrack(const std::vector<std::pair<ViewId, Feature> >& track);

  // Removes the track from the reconstruction including the corresponding
  // features that are present in the view that observe it.
  bool RemoveTrack(const TrackId track_id);
  int NumTracks() const;

  // Returns the Track or a nullptr if the track does not exist.
  const class Track* Track(const TrackId track_id) const;
  class Track* MutableTrack(const TrackId track_id);

  // Return all TrackIds in the reconstruction.
  std::vector<TrackId> TrackIds() const;

  // Normalizes the reconstruction such that the "center" of the reconstruction
  // is moved to the origin and the reconstruction is scaled such that the
  // median distance of 3D points from the origin is 100.0. This does not affect
  // the reprojection error. A rotation is applied such that the x-z plane is
  // set to the dominating plane of the cameras.
  //
  // NOTE: This implementation is inspired by the BAL problem normalization in
  // Ceres Solver.
  void Normalize();

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar) {  // NOLINT
    ar(next_track_id_,
       next_view_id_,
       view_name_to_id_,
       views_,
       tracks_);
  }

  TrackId next_track_id_;
  ViewId next_view_id_;
  CameraIntrinsicsGroupId next_camera_intrinsics_group_id_;

  std::unordered_map<std::string, ViewId> view_name_to_id_;
  std::unordered_map<ViewId, class View> views_;
  std::unordered_map<TrackId, class Track> tracks_;

  std::unordered_map<ViewId, CameraIntrinsicsGroupId>
      view_id_to_camera_intrinsics_group_id_;
  std::unordered_map<CameraIntrinsicsGroupId, std::unordered_set<ViewId> >
      camera_intrinsics_groups_;
};

}  // namespace theia

#endif  // THEIA_SFM_RECONSTRUCTION_H_
