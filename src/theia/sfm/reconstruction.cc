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

#include "theia/sfm/reconstruction.h"

#include <Eigen/Core>
#include <Eigen/SVD>
#include <glog/logging.h>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "theia/util/map_util.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/pose/util.h"
#include "theia/sfm/track.h"
#include "theia/sfm/transformation/transform_reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {

namespace {

// Return whether the track contains the same view twice.
bool DuplicateViewsExistInTrack(
    const std::vector<std::pair<ViewId, Feature> >& track) {
  std::vector<ViewId> view_ids;
  view_ids.reserve(track.size());
  for (const auto& feature : track) {
    view_ids.push_back(feature.first);
  }
  std::sort(view_ids.begin(), view_ids.end());
  return (std::adjacent_find(view_ids.begin(), view_ids.end()) !=
          view_ids.end());
}

double Median(std::vector<double>* data) {
  int n = data->size();
  std::vector<double>::iterator mid_point = data->begin() + n / 2;
  std::nth_element(data->begin(), mid_point, data->end());
  return *mid_point;
}

}  // namespace

Reconstruction::Reconstruction()
    : next_track_id_(0),
      next_view_id_(0),
      next_camera_intrinsics_group_id_(0) {}

Reconstruction::~Reconstruction() {}

ViewId Reconstruction::ViewIdFromName(const std::string& view_name) const {
  return FindWithDefault(view_name_to_id_, view_name, kInvalidViewId);
}

ViewId Reconstruction::AddView(const std::string& view_name) {
  const ViewId view_id = AddView(view_name, next_camera_intrinsics_group_id_);
  ++next_camera_intrinsics_group_id_;
  return view_id;
}

ViewId Reconstruction::AddView(const std::string& view_name,
                               const CameraIntrinsicsGroupId group_id) {
  if (ContainsKey(view_name_to_id_, view_name)) {
    LOG(WARNING) << "Could not add view with the name " << view_name
               << " because that name already exists in the reconstruction.";
    return kInvalidViewId;
  }

  if (view_name.empty()) {
    LOG(WARNING)
        << "View name was empty. Could not add view to reconstruction.";
    return kInvalidViewId;
  }

  class View new_view(view_name);

  // If the camera intrinsics group already exists, set the internal intrinsics
  // of each camera to point to the same underlying intrinsics.
  if (camera_intrinsics_groups_[group_id].size() > 0) {
    const ViewId view_id_in_intrinsics_group =
        *camera_intrinsics_groups_[group_id].begin();
    const Camera& intrinsics_group_camera =
        FindOrDie(views_, view_id_in_intrinsics_group).Camera();
    LOG(INFO) << "Attempting to set the intrinsics of view " << next_view_id_
              << " and " << view_id_in_intrinsics_group
              << " to shared intrinsics.";

    // Set the shared_ptr objects to point to the same place so that the
    // intrinsics are truly shared.
    new_view.MutableCamera()->MutableCameraIntrinsics() =
        intrinsics_group_camera.CameraIntrinsics();
    LOG(INFO)
        << "New view location: "
        << new_view.MutableCamera()->MutableCameraIntrinsics()->parameters()
        << " vs old location: "
        << intrinsics_group_camera.CameraIntrinsics()->parameters();
  }

  // Add the view to the reconstruction.
  views_.emplace(next_view_id_, new_view);
  view_name_to_id_.emplace(view_name, next_view_id_);

  // Add this view to the camera intrinsics group, and vice versa.
  view_id_to_camera_intrinsics_group_id_.emplace(next_view_id_, group_id);
  camera_intrinsics_groups_[group_id].emplace(next_view_id_);

  ++next_view_id_;
  return next_view_id_ - 1;
}

bool Reconstruction::RemoveView(const ViewId view_id) {
  class View* view = FindOrNull(views_, view_id);
  if (view == nullptr) {
    LOG(WARNING)
        << "Could not remove the view from the reconstruction because the view "
           "does not exist.";
    return false;
  }

  const auto& track_ids = view->TrackIds();
  for (const TrackId track_id : track_ids) {
    class Track* track = MutableTrack(track_id);
    if (track == nullptr) {
      LOG(WARNING) << "Could not remove the view from the track because the "
                    "track does not exist";
      return false;
    }

    if (!track->RemoveView(view_id)) {
      LOG(WARNING) << "Could not remove to view from the track";
      return false;
    }

    // Remove the track if it no longer contains any views.
    if (track->NumViews() == 0) {
      RemoveTrack(track_id);
    }
  }

  // Remove the view name.
  const std::string& view_name = view->Name();
  view_name_to_id_.erase(view_name);

  // Remove the view from the camera intrinsics groups.
  const CameraIntrinsicsGroupId group_id =
      CameraIntrinsicsGroupIdFromViewId(view_id);
  view_id_to_camera_intrinsics_group_id_.erase(view_id);
  std::unordered_set<ViewId>& camera_intrinsics_group =
      FindOrDie(camera_intrinsics_groups_, group_id);
  camera_intrinsics_group.erase(view_id);
  if (camera_intrinsics_group.size() == 0) {
    camera_intrinsics_groups_.erase(group_id);
  }

  // Remove the view.
  views_.erase(view_id);
  return true;
}

int Reconstruction::NumViews() const {
  return views_.size();
}

const class View* Reconstruction::View(const ViewId view_id) const {
  return FindOrNull(views_, view_id);
}

class View* Reconstruction::MutableView(const ViewId view_id) {
  return FindOrNull(views_, view_id);
}

std::vector<ViewId> Reconstruction::ViewIds() const {
  std::vector<ViewId> view_ids;
  view_ids.reserve(views_.size());
  for (const auto& view : views_) {
    view_ids.push_back(view.first);
  }
  return view_ids;
}

  // Get the camera intrinsics group id for the view id.
CameraIntrinsicsGroupId Reconstruction::CameraIntrinsicsGroupIdFromViewId(
    const ViewId view_id) const {
  return FindWithDefault(view_id_to_camera_intrinsics_group_id_,
                         view_id,
                         kInvalidCameraIntrinsicsGroupId);
}

// Return all view ids with the given camera intrinsics group id. If an
// invalid or non-existant group is chosen then an empty set will be returned.
std::unordered_set<ViewId> Reconstruction::GetViewsInCameraIntrinsicGroup(
    const CameraIntrinsicsGroupId group_id) const {
  return FindWithDefault(camera_intrinsics_groups_,
                         group_id,
                         std::unordered_set<ViewId>());
}

std::unordered_set<CameraIntrinsicsGroupId>
Reconstruction::CameraIntrinsicsGroupIds() const {
  std::unordered_set<CameraIntrinsicsGroupId> group_ids;
  group_ids.reserve(camera_intrinsics_groups_.size());
  for (const auto& groups : camera_intrinsics_groups_) {
    group_ids.emplace(groups.first);
  }
  return group_ids;
}

int Reconstruction::NumCameraIntrinsicGroups() const {
  return camera_intrinsics_groups_.size();
}

TrackId Reconstruction::AddTrack() {
  const TrackId new_track_id = next_track_id_;
  CHECK(!ContainsKey(tracks_, new_track_id))
      << "The reconstruction already contains a track with id: "
      << new_track_id;

  class Track new_track;
  tracks_.emplace(new_track_id, new_track);
  ++next_track_id_;
  return new_track_id;
}

bool Reconstruction::AddObservation(const ViewId view_id,
                                    const TrackId track_id,
                                    const Feature& feature) {
  CHECK(ContainsKey(views_, view_id))
      << "View does not exist. AddObservation may only be used to add "
         "observations to an existing view.";
  CHECK(ContainsKey(tracks_, track_id))
      << "Track does not exist. AddObservation may only be used to add "
         "observations to an existing track.";

  class View* view = FindOrNull(views_, view_id);
  class Track* track = FindOrNull(tracks_, track_id);
  if (view->GetFeature(track_id) != nullptr) {
    LOG(WARNING)
        << "Cannot add a new observation of track " << track_id
        << " because the view already contains an observation of the track.";
    return false;
  }

  const std::unordered_set<ViewId>& views_observing_track = track->ViewIds();
  if (ContainsKey(views_observing_track, view_id)) {
    LOG(WARNING)
        << "Cannot add a new observation of track " << track_id
        << " because the track is already observed by view " << view_id;
    return false;
  }

  view->AddFeature(track_id, feature);
  track->AddView(view_id);
  return true;
}

TrackId Reconstruction::AddTrack(
    const std::vector<std::pair<ViewId, Feature> >& track) {
  if (track.size() < 2) {
    LOG(WARNING) << "Tracks must have at least 2 observations (" << track.size()
                 << " were given). Cannot add track to the reconstruction";
    return kInvalidTrackId;
  }

  if (DuplicateViewsExistInTrack(track)) {
    LOG(WARNING) << "Cannot add a track that contains the same view twice to "
                    "the reconstruction.";
    return kInvalidTrackId;
  }

  const TrackId new_track_id = next_track_id_;
  CHECK(!ContainsKey(tracks_, new_track_id))
      << "The reconstruction already contains a track with id: "
      << new_track_id;

  class Track new_track;
  for (const auto& observation : track) {
    // Make sure the view exists in the model.
    CHECK(ContainsKey(views_, observation.first))
        << "Cannot add a track with containing an observation in view id "
        << observation.first << " because the view does not exist.";

    // Add view to track.
    new_track.AddView(observation.first);

    // Add track to view.
    class View* view = MutableView(observation.first);
    view->AddFeature(new_track_id, observation.second);
  }

  tracks_.emplace(new_track_id, new_track);
  ++next_track_id_;
  return new_track_id;
}

bool Reconstruction::RemoveTrack(const TrackId track_id) {
  class Track* track = FindOrNull(tracks_, track_id);
  if (track == nullptr) {
    LOG(WARNING) << "Cannot remove a track that does not exist";
    return false;
  }

  // Remove track from views.
  for (const ViewId view_id : track->ViewIds()) {
    class View* view = FindOrNull(views_, view_id);
    if (view == nullptr) {
      LOG(WARNING) << "Could not remove a track from the view because the view "
                    "does not exist";
      return false;
    }

    if (!view->RemoveFeature(track_id)) {
      LOG(WARNING) << "Could not remove the track from the view because the "
                    "track is not observed by the view.";
      return false;
    }
  }

  // Delete from the reconstruction.
  tracks_.erase(track_id);
  return true;
}

int Reconstruction::NumTracks() const {
  return tracks_.size();
}

const class Track* Reconstruction::Track(const TrackId track_id) const {
  return FindOrNull(tracks_, track_id);
}

class Track* Reconstruction::MutableTrack(const TrackId track_id) {
  return FindOrNull(tracks_, track_id);
}

std::vector<TrackId> Reconstruction::TrackIds() const {
  std::vector<TrackId> track_ids;
  track_ids.reserve(tracks_.size());
  for (const auto& track : tracks_) {
    track_ids.push_back(track.first);
  }
  return track_ids;
}

void Reconstruction::Normalize() {
  // First normalize the position so that the marginal median of the camera
  // positions is at the origin.
  std::vector<std::vector<double> > camera_positions(3);
  Eigen::Vector3d median_camera_position;
  const auto& view_ids = ViewIds();
  for (const ViewId view_id : view_ids) {
    const class View* view = View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }

    const Eigen::Vector3d point = view->Camera().GetPosition();
    camera_positions[0].push_back(point[0]);
    camera_positions[1].push_back(point[1]);
    camera_positions[2].push_back(point[2]);
  }
  median_camera_position(0) = Median(&camera_positions[0]);
  median_camera_position(1) = Median(&camera_positions[1]);
  median_camera_position(2) = Median(&camera_positions[2]);
  TransformReconstruction(Eigen::Matrix3d::Identity(), -median_camera_position,
                          1.0, this);

  // Get the estimated track ids.
  const auto& temp_track_ids = TrackIds();
  std::vector<TrackId> track_ids;
  track_ids.reserve(temp_track_ids.size());
  for (const TrackId track_id : temp_track_ids) {
    const class Track* track = Track(track_id);
    if (track == nullptr || !track->IsEstimated()) {
      continue;
    }
    track_ids.emplace_back(track_id);
  }

  if (track_ids.size() == 0) {
    return;
  }

  // Compute the marginal median of the 3D points.
  std::vector<std::vector<double> > points(3);
  Eigen::Vector3d median;
  for (const TrackId track_id : track_ids) {
    const class Track* track = Track(track_id);
    const Eigen::Vector3d point = track->Point().hnormalized();
    points[0].push_back(point[0]);
    points[1].push_back(point[1]);
    points[2].push_back(point[2]);
  }
  median(0) = Median(&points[0]);
  median(1) = Median(&points[1]);
  median(2) = Median(&points[2]);

  // Find the median absolute deviation of the points from the median.
  std::vector<double> distance_to_median;
  distance_to_median.reserve(track_ids.size());
  for (const TrackId track_id : track_ids) {
    const class Track* track = Track(track_id);
    const Eigen::Vector3d point = track->Point().hnormalized();
    distance_to_median.emplace_back((point - median).lpNorm<1>());
  }
  // This will scale the reconstruction so that the median absolute deviation of
  // the points is 100.
  const double scale = 100.0 / Median(&distance_to_median);
  TransformReconstruction(Eigen::Matrix3d::Identity(),
                          Eigen::Vector3d::Zero(),
                          scale,
                          this);

  // Most images are taken relatively upright with the x-direction of the image
  // parallel to the ground plane. We can solve for the transformation that
  // tries to best align the x-directions to the ground plane by finding the
  // null vector of the covariance matrix of per-camera x-directions.
  Eigen::Matrix3d correlation;
  correlation.setZero();
  for (const ViewId view_id : view_ids) {
    const class View* view = View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }
    const Camera& camera = View(view_id)->Camera();
    const Eigen::Vector3d x =
        camera.GetOrientationAsRotationMatrix().transpose() *
        Eigen::Vector3d(1.0, 0.0, 0.0);
    correlation += x * x.transpose();
  }

  // The up-direction is computed as the null vector of the covariance matrix.
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(correlation, Eigen::ComputeFullV);
  const Eigen::Vector3d plane_normal = svd.matrixV().rightCols<1>();

  // We want the coordinate system to be such that the cameras lie on the x-z
  // plane with the y vector pointing up. Thus, the plane normal should be equal
  // to the positive y-direction.
  Eigen::Matrix3d rotation = Eigen::Quaterniond::FromTwoVectors(
      plane_normal, Eigen::Vector3d(0, 1, 0)).toRotationMatrix();
  TransformReconstruction(rotation, Eigen::Vector3d::Zero(), 1.0, this);
}

}  // namespace theia
