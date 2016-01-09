// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include <Eigen/Core>
#include <glog/logging.h>
#include <string>
#include <utility>
#include <vector>

#include "theia/io/import_nvm_file.h"
#include "theia/sfm/reconstruction.h"
#include "theia/util/filesystem.h"
#include "theia/util/map_util.h"

#include "visual_sfm/util.h"
#include "visual_sfm/DataInterface.h"

namespace theia {

namespace {
struct VSFMReconstruction {
  std::vector<CameraT> camera_data;
  std::vector<Point3D> point_data;
  std::vector<Point2D> measurements;
  std::vector<int> point_index;
  std::vector<int> camera_index;
  std::vector<string> view_names;
  std::vector<int> point_colors;
};
}  // namespace

bool ImportNVMFile(const std::string& nvm_filepath,
                   Reconstruction* reconstruction) {
  CHECK_GT(nvm_filepath.length(), 0);
  CHECK_NOTNULL(reconstruction);

  // Read VSFM file.
  VSFMReconstruction vsfm_reconstruction;
  CHECK(LoadModelFile(nvm_filepath.c_str(),
                      vsfm_reconstruction.camera_data,
                      vsfm_reconstruction.point_data,
                      vsfm_reconstruction.measurements,
                      vsfm_reconstruction.point_index,
                      vsfm_reconstruction.camera_index,
                      vsfm_reconstruction.view_names,
                      vsfm_reconstruction.point_colors));

  // Check the NVM reconstruction for sanity.
  CHECK_EQ(vsfm_reconstruction.view_names.size(),
           vsfm_reconstruction.camera_data.size());
  CHECK_EQ(vsfm_reconstruction.measurements.size(),
           vsfm_reconstruction.point_index.size());
  CHECK_EQ(vsfm_reconstruction.point_index.size(),
           vsfm_reconstruction.camera_index.size());

  // Add all cameras to the reconstruction.
  for (int i = 0; i < vsfm_reconstruction.camera_data.size(); i++) {
    std::string view_name;
    GetFilenameFromFilepath(vsfm_reconstruction.view_names[i],
                                   true,
                                   &view_name);
    // Add the view to the reconstruction.
    LOG(INFO) << "Adding view " << view_name << " to the reconstruction.";
    const ViewId view_id = reconstruction->AddView(view_name);
    CHECK_NE(view_id, kInvalidViewId);
    View* view = reconstruction->MutableView(view_id);
    view->SetEstimated(true);

    // Set the camera intrinsic parameters.
    const CameraT& vsfm_camera = vsfm_reconstruction.camera_data[i];
    Camera* camera = view->MutableCamera();
    camera->SetFocalLength(vsfm_camera.GetFocalLength());

    // Set the camera extrinsic parameters. The rotation matrix is retreived
    // from visual sfm in row-major order, hence the transpose.
    Eigen::Matrix3f rotation;
    vsfm_camera.GetMatrixRotation(rotation.data());
    camera->SetOrientationFromRotationMatrix(
        rotation.transpose().cast<double>());
    Eigen::Vector3f position;
    vsfm_camera.GetCameraCenter(position.data());
    camera->SetPosition(position.cast<double>());
  }

  // Create the track correspondences.
  std::unordered_map<int, std::vector<std::pair<ViewId, Feature>>>
      tracks;
  for (int i = 0; i < vsfm_reconstruction.measurements.size(); i++) {
    const Feature feature(vsfm_reconstruction.measurements[i].x,
                          vsfm_reconstruction.measurements[i].y);
    const int img_id = vsfm_reconstruction.camera_index[i];
    tracks[vsfm_reconstruction.point_index[i]].emplace_back(img_id, feature);
  }

  // Add all tracks to the reconstruction and set the 3d position.
  for (int i = 0; i < vsfm_reconstruction.point_data.size(); i++) {
    Eigen::Vector4f point(vsfm_reconstruction.point_data[i].xyz[0],
                          vsfm_reconstruction.point_data[i].xyz[1],
                          vsfm_reconstruction.point_data[i].xyz[2],
                          1.0);
    Eigen::Vector3i color(vsfm_reconstruction.point_colors[3 * i],
                          vsfm_reconstruction.point_colors[3 * i + 1],
                          vsfm_reconstruction.point_colors[3 * i + 2]);
    const auto& features = FindOrDie(tracks, i);
    const TrackId track_id = reconstruction->AddTrack(features);
    Track* track = reconstruction->MutableTrack(track_id);
    track->SetEstimated(true);
    *track->MutablePoint() = point.cast<double>();
    *track->MutableColor() = color.cast<uint8_t>();
  }

  return true;
}

}  // namespace theia
