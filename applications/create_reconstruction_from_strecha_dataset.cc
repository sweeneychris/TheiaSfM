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
//
// This program is a useful tool for converting camera files from the Strecha
// Multi-View Stereo dataset into Theia reconstructions. The output
// reconstruction may serve as ground truth for experiments that utilize the
// Theia library. The datasets may be found at the link below, and the relevant
// files for this program are the "camera" files.
//
// http://cvlabwww.epfl.ch/data/multiview/denseMVS.html

#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>

#include <fstream>  // NOLINT
#include <string>

DEFINE_string(input_dataset_directory, "",
              "Directory containing camera files from the Strecha MVS dataset. "
              "Do not include a trailing slash.");
DEFINE_string(output_reconstruction, "",
              "Filepath of the output Theia reconstruction file generated.");

using theia::Camera;
using theia::PinholeCameraModel;
using theia::Reconstruction;
using theia::ViewId;

void AddCameraToReconstruction(const std::string& camera_filepath,
                               Reconstruction* reconstruction) {
  std::ifstream ifs(camera_filepath, std::ios::in);
  if (!ifs.is_open()) {
    LOG(FATAL) << "Cannot read the camera file from: " << camera_filepath;
  }
  LOG(INFO) << "Reading camera parameters from file: " << camera_filepath;

  std::string image_name;
  CHECK(theia::GetFilenameFromFilepath(camera_filepath, true, &image_name));

  // Remove the ".camera" from the filename so we are just left with the image
  // name.
  const std::size_t dot_camera_pos = image_name.find(".camera");
  image_name = image_name.substr(0, dot_camera_pos);

  const ViewId view_id = reconstruction->AddView(image_name);
  CHECK_NE(view_id, theia::kInvalidViewId)
      << "The image " << image_name
      << " could not be added to the reconstruction.";

  Camera* camera = reconstruction->MutableView(view_id)->MutableCamera();

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
  theia::CalibrationMatrixToIntrinsics(
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

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  // Ensure the input parameters are well-formed.
  CHECK_GT(FLAGS_input_dataset_directory.length(), 0);
  CHECK_GT(FLAGS_output_reconstruction.length(), 0);

  // Get the filepaths of the Strecha MVS dataset camera files.
  const std::string camera_wildcard =
      FLAGS_input_dataset_directory + "/*.camera";
  std::vector<std::string> camera_files;
  CHECK(theia::GetFilepathsFromWildcard(camera_wildcard, &camera_files))
      << "Could not find cameras that matched the filepath: " << camera_wildcard
      << ". NOTE that the ~ filepath is not supported.";

  CHECK_GT(camera_files.size(), 0)
      << "No cameras found in: " << camera_wildcard;

  // Add each camera in the Strecha MVS dataset to the reconstruciton.
  Reconstruction reconstruction;
  for (const std::string camera_file : camera_files) {
    AddCameraToReconstruction(camera_file, &reconstruction);
  }

  CHECK(theia::WriteReconstruction(reconstruction, FLAGS_output_reconstruction))
      << "Could not write the reconstruction file to: "
      << FLAGS_output_reconstruction;

  return 0;
}
