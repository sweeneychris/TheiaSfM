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
#include "theia/sfm/estimators/estimate_dominant_plane_from_points.h"
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
  views_.emplace(next_view_id_, new_view);
  view_name_to_id_.emplace(view_name, next_view_id_);

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

  std::vector<ViewId> views_to_add;
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

  // Apply position transformation.
  TransformReconstruction(Eigen::Matrix3d::Identity(), -median, 1.0, this);

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

  // Compute a rotation such that the x-z plane is aligned to the dominating
  // plane of the cameras.
  std::vector<Eigen::Vector3d> cameras;
  const auto& view_ids = this->ViewIds();
  for (const ViewId view_id : view_ids) {
    const class View* view = View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }
    cameras.emplace_back(view->Camera().GetPosition());
  }

  // Robustly estimate the dominant plane from the cameras. This will correspond
  // to a plan that is parallel to the ground plane for the majority of
  // reconstructions. We start with a small threshold and gradually increase it
  // until at an inlier set of at least 50% is found.
  RansacParameters ransac_params;
  ransac_params.max_iterations = 1000;
  ransac_params.error_thresh = 0.01;
  Plane plane;
  RansacSummary unused_summary;
  Eigen::Matrix3d rotation_for_dominant_plane = Eigen::Matrix3d::Identity();
  if (EstimateDominantPlaneFromPoints(ransac_params,
                                      RansacType::LMED,
                                      cameras,
                                      &plane,
                                      &unused_summary)) {
    // Set the rotation such that the plane normal points in the upward
    // direction. Choose the sign of the normal that will minimize the rotation
    // (this hopes to prevent having a rotation that flips the scene upside
    // down).
    const Eigen::Quaterniond rotation_quat1 =
        Eigen::Quaterniond::FromTwoVectors(plane.unit_normal,
                                           Eigen::Vector3d(0, 1.0, 0));
    const Eigen::Quaterniond rotation_quat2 =
        Eigen::Quaterniond::FromTwoVectors(-plane.unit_normal,
                                           Eigen::Vector3d(0, 1.0, 0));
    const Eigen::AngleAxisd rotation1_aa(rotation_quat1);
    const Eigen::AngleAxisd rotation2_aa(rotation_quat2);

    if (rotation1_aa.angle() < rotation2_aa.angle()) {
      rotation_for_dominant_plane = rotation1_aa.toRotationMatrix();
    } else {
      rotation_for_dominant_plane = rotation2_aa.toRotationMatrix();
    }
  }

  TransformReconstruction(rotation_for_dominant_plane,
                          Eigen::Vector3d::Zero(),
                          1.0,
                          this);
}

}  // namespace theia
