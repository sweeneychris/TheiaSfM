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
// Author: Aleksander Holynski (holynski@cs.washington.edu)

#include "theia/io/write_colmap_files.h"

#include <Eigen/Core>
#include <fstream>  // NOLINT
#include <glog/logging.h>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/camera/pinhole_camera_model.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"
#include "theia/util/map_util.h"

namespace theia {
namespace {

bool WriteCamerasFile(const Reconstruction& reconstruction,
                      const std::string& cameras_file) {
  std::ofstream ofs_cameras(cameras_file);
  if (!ofs_cameras.is_open()) {
    LOG(ERROR) << "Cannot open the file: " << cameras_file << " for writing.";
    return false;
  }

  const auto& group_ids = reconstruction.CameraIntrinsicsGroupIds();

  for (auto group_id : group_ids) {
    const View* view = reconstruction.View(
        *reconstruction.GetViewsInCameraIntrinsicGroup(group_id).begin());
    const Camera& camera = view->Camera();
    if (camera.GetCameraIntrinsicsModelType() !=
        CameraIntrinsicsModelType::PINHOLE) {
      LOG(ERROR) << "Could not add camera " << view->Name()
                 << " to the COLMAP output file because COLMAP doesn't support "
                    "FOV or fisheye cameras.";
      return false;
    }

    const CameraIntrinsicsModel& intrinsics = *camera.CameraIntrinsics();
    ofs_cameras
        << group_id << " RADIAL " << camera.ImageWidth() << " "
        << camera.ImageHeight() << " " << camera.FocalLength() << " "
        << camera.PrincipalPointX() << " " << camera.PrincipalPointY() << " "
        << intrinsics.GetParameter(PinholeCameraModel::RADIAL_DISTORTION_1)
        << " "
        << intrinsics.GetParameter(PinholeCameraModel::RADIAL_DISTORTION_1)
        << std::endl;
  }
  ofs_cameras.close();

  return true;
}

bool WriteImagesFile(const Reconstruction& reconstruction,
                     const std::string& images_file) {
  std::ofstream ofs_images(images_file);

  if (!ofs_images.is_open()) {
    LOG(ERROR) << "Cannot open the file: " << images_file << " for writing.";
    return false;
  }

  for (auto view_id : reconstruction.ViewIds()) {
    const View* view = reconstruction.View(view_id);
    const Camera& camera = view->Camera();
    const Eigen::Vector3d translation =
        -camera.GetOrientationAsRotationMatrix() * camera.GetPosition();
    const Eigen::Quaterniond orientation(
        camera.GetOrientationAsRotationMatrix());
    ofs_images << view_id << " " << orientation.w() << " " << orientation.x()
               << " " << orientation.y() << " " << orientation.z() << " "
               << translation.x() << " " << translation.y() << " "
               << translation.z() << " "
               << reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id)
               << " " << view->Name() << std::endl;
    const auto& track_ids = view->TrackIds();
    for (auto track_id : track_ids) {
      auto feature = view->GetFeature(track_id);
      ofs_images << feature->x() << " " << feature->y() << " " << track_id
                 << " ";
    }
    ofs_images << std::endl;
  }
  ofs_images.close();
  return true;
}

bool WritePointsFile(const Reconstruction& reconstruction,
                     const std::string& points_file) {
  std::ofstream ofs_points(points_file);
  if (!ofs_points.is_open()) {
    LOG(ERROR) << "Cannot open the file: " << points_file << " for writing.";
    return false;
  }

  for (auto track_id : reconstruction.TrackIds()) {
    const Track* track = reconstruction.Track(track_id);
    const Eigen::Vector3d point = track->Point().hnormalized();
    const auto& color = track->Color();
    ofs_points << track_id << " " << point.x() << " " << point.y() << " "
               << point.z() << " " << static_cast<int>(color[0]) << " "
               << static_cast<int>(color[1]) << " "
               << static_cast<int>(color[2]) << " " << 0.0 << " ";
    for (auto view_id : track->ViewIds()) {
      const View* view = reconstruction.View(view_id);
      const auto& track_ids = view->TrackIds();
      const int point_index =
          std::find(track_ids.begin(), track_ids.end(), track_id) -
          track_ids.begin();
      ofs_points << view_id << " " << point_index << " ";
    }
    ofs_points << std::endl;
  }

  ofs_points.close();
  return true;
}

}  // namespace

bool WriteColmapFiles(const Reconstruction& reconstruction,
                      const std::string& output_directory) {
  Reconstruction estimated_reconstruction;
  CreateEstimatedSubreconstruction(reconstruction, &estimated_reconstruction);

  const std::string& cameras_file = output_directory + "/cameras.txt";
  const std::string& images_file = output_directory + "/images.txt";
  const std::string& points_file = output_directory + "/points3D.txt";

  if (!WriteCamerasFile(estimated_reconstruction, cameras_file)) {
    return false;
  }
  if (!WriteImagesFile(estimated_reconstruction, images_file)) {
    return false;
  }
  if (!WritePointsFile(estimated_reconstruction, points_file)) {
    return false;
  }
  return true;
}

}  // namespace theia
