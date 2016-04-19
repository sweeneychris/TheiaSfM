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

#include "theia/sfm/view.h"

#include <string>
#include <unordered_set>
#include <unordered_map>

#include "theia/util/map_util.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/types.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {

View::View()
    : name_(""), is_estimated_(false),
      camera_intrinsics_group_id_(kInvalidCameraIntrinsicsGroupId) {}

View::View(const std::string &name)
    : name_(name), is_estimated_(false),
      camera_intrinsics_group_id_(kInvalidCameraIntrinsicsGroupId) {}

const std::string& View::Name() const {
  return name_;
}

void View::SetEstimated(bool is_estimated) {
  is_estimated_ = is_estimated;
}

bool View::IsEstimated() const {
  return is_estimated_;
}

const class Camera& View::Camera() const {
  return camera_;
}

class Camera* View::MutableCamera() {
  return &camera_;
}

CameraIntrinsicsGroupId View::GetCameraIntrinsicsGroupId() const {
  return camera_intrinsics_group_id_;
}

void View::SetCameraIntrinsicsGroupId(const CameraIntrinsicsGroupId group_id) {
  // TODO(stoyanovd): is this check needed here?
  CHECK(group_id != kInvalidCameraIntrinsicsGroupId)
    << "Could not set " << "camera_intrinsics_group_id to invalid value.";
  camera_intrinsics_group_id_ = group_id;
}

const struct CameraIntrinsicsPrior& View::CameraIntrinsicsPrior() const {
  return camera_intrinsics_prior_;
}

struct CameraIntrinsicsPrior* View::MutableCameraIntrinsicsPrior() {
  return &camera_intrinsics_prior_;
}

int View::NumFeatures() const {
  return features_.size();
}

std::vector<TrackId> View::TrackIds() const {
  std::vector<TrackId> track_ids;
  track_ids.reserve(features_.size());
  for (const auto& track : features_) {
    track_ids.emplace_back(track.first);
  }
  return track_ids;
}

const Feature* View::GetFeature(const TrackId track_id) const {
  return FindOrNull(features_, track_id);
}

void View::AddFeature(const TrackId track_id, const Feature& feature) {
  features_[track_id] = feature;
}

bool View::RemoveFeature(const TrackId track_id) {
  return features_.erase(track_id) > 0;
}

}  // namespace theia
