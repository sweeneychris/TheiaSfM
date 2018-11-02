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

#include "theia/io/read_strecha_dataset.h"

#include <fstream>
#include <iostream>
#include <string>

#include "theia/sfm/camera/pinhole_camera_model.h"
#include "theia/sfm/camera/projection_matrix_utils.h"
#include "theia/sfm/reconstruction.h"
#include "theia/util/filesystem.h"
#include "theia/util/string.h"

namespace theia {

void AddCameraToReconstruction(const std::string& camera_filepath,
                               Reconstruction* reconstruction) {
  static const int kCameraIntrinsicsGroup = 0;

  std::ifstream ifs(camera_filepath, std::ios::in);
  if (!ifs.is_open()) {
    LOG(FATAL) << "Cannot read the camera file from: " << camera_filepath;
  }
  LOG(INFO) << "Reading camera parameters from file: " << camera_filepath;

  std::string image_name;
  CHECK(GetFilenameFromFilepath(camera_filepath, true, &image_name));

  // Remove the ".camera" from the filename so we are just left with the image
  // name.
  const std::size_t dot_camera_pos = image_name.find(".camera");
  image_name = image_name.substr(0, dot_camera_pos);

  // Add the view and ensure that the camera intrinsics are all the same.
  const ViewId view_id =
      reconstruction->AddView(image_name, kCameraIntrinsicsGroup);
  CHECK_NE(view_id, kInvalidViewId)
      << "The image " << image_name
      << " could not be added to the reconstruction.";

  Camera* camera = reconstruction->MutableView(view_id)->MutableCamera();
  camera->SetCameraIntrinsicsModelType(CameraIntrinsicsModelType::PINHOLE);

  // The first three rows are the calibration matrix.
  Eigen::Matrix3d calibration_matrix;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ifs >> calibration_matrix(i, j);
    }
  }
  // Extract the camera intrinsics from the calibration matrix.
  double* camera_intrinsics =
      camera->MutableCameraIntrinsics()->mutable_parameters();
  CalibrationMatrixToIntrinsics(
      calibration_matrix,
      camera_intrinsics + PinholeCameraModel::FOCAL_LENGTH,
      camera_intrinsics + PinholeCameraModel::SKEW,
      camera_intrinsics + PinholeCameraModel::ASPECT_RATIO,
      camera_intrinsics + PinholeCameraModel::PRINCIPAL_POINT_X,
      camera_intrinsics + PinholeCameraModel::PRINCIPAL_POINT_Y);

  // The next line is simply 0 0 0.
  int unused_zero;
  ifs >> unused_zero >> unused_zero >> unused_zero;

  // Then comes the rotation matrix.
  Eigen::Matrix3d rotation_matrix;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ifs >> rotation_matrix(j, i);
    }
  }
  camera->SetOrientationFromRotationMatrix(rotation_matrix);

  // Next, the camera position.
  Eigen::Vector3d position;
  ifs >> position(0) >> position(1) >> position(2);
  camera->SetPosition(position);

  // Image size (width, height).
  int image_width, image_height;
  ifs >> image_width >> image_height;
  camera->SetImageSize(image_width, image_height);

  // Set the camera to estimated, otherwise it will not be saved when writing
  // the reconstruction.
  reconstruction->MutableView(view_id)->SetEstimated(true);
}

// Loads all camera intrinsics and extrinsics information from the Strecha MVS
// dataset. A Reconstruction object is populated with the camera information.
bool ReadStrechaDataset(const std::string& dataset_directory,
                        Reconstruction* reconstruction) {
  // Get the filepaths of the Strecha MVS dataset camera files.
  std::string camera_wildcard = dataset_directory;
  AppendTrailingSlashIfNeeded(&camera_wildcard);
  camera_wildcard = camera_wildcard + "/*.camera";

  std::vector<std::string> camera_files;
  CHECK(GetFilepathsFromWildcard(camera_wildcard, &camera_files))
      << "Could not find cameras that matched the filepath: " << camera_wildcard
      << ". NOTE that the ~ filepath is not supported.";

  CHECK_GT(camera_files.size(), 0)
      << "No cameras found in: " << camera_wildcard;

  // Add each camera in the Strecha MVS dataset to the reconstruciton.
  for (const std::string camera_file : camera_files) {
    AddCameraToReconstruction(camera_file, reconstruction);
  }
  return true;
}

}  // namespace theia
