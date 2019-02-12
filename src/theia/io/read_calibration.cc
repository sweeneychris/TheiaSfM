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

#include "theia/io/read_calibration.h"

#include <stdio.h>
#include <glog/logging.h>
#include <string>
#include <unordered_map>

#include <rapidjson/document.h>

#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/camera/camera_intrinsics_model_type.h"

namespace theia {
namespace {

static const char* kPinholeType = "PINHOLE";
static const char* kPriorsEntry = "priors";
static const char* kCameraIntrinsicsPrior = "CameraIntrinsicsPrior";
static const char* kImageName = "image_name";
static const char* kCameraType = "camera_intrinsics_type";

// Pinhole camera parameters.
static const char* kFocalLength = "focal_length";
static const char* kImageWidth = "width";
static const char* kImageHeight = "height";
static const char* kPrincipalPointX = "principal_point_x";
static const char* kPrincipalPointY = "principal_point_y";
static const char* kAspectRatio = "aspect_ratio";
static const char* kSkew = "skew";
static const char* kRadialDistortionCoeffs = "radial_distortion_coeffs";

bool ExtractPinholeCamera(const rapidjson::Value& entry,
                          CameraIntrinsicsPrior* prior) {
  // Get the focal length.
  if (entry.HasMember(kFocalLength)) {
    prior->focal_length.is_set = true;
    prior->focal_length.value[0] = entry[kFocalLength].GetDouble();
  }

  // Get the principal points.
  if (entry.HasMember(kPrincipalPointX) && entry.HasMember(kPrincipalPointY)) {
    prior->principal_point.is_set = true;
    prior->principal_point.value[0] = entry[kPrincipalPointX].GetDouble();
    prior->principal_point.value[1] = entry[kPrincipalPointY].GetDouble();
    prior->image_width = static_cast<int>(2 * prior->principal_point.value[0]);
    prior->image_height = static_cast<int>(2 * prior->principal_point.value[1]);
  }

  // Get width.
  if (entry.HasMember(kImageWidth)) {
    prior->image_width = entry[kImageWidth].GetInt();
  }

  // Get height.
  if (entry.HasMember(kImageHeight)) {
    prior->image_height = entry[kImageHeight].GetInt();
  }

  // Get aspect ratio.
  if (entry.HasMember(kAspectRatio)) {
    prior->aspect_ratio.is_set = true;
    prior->aspect_ratio.value[0] = entry[kAspectRatio].GetDouble();
  }

  // Get skew.
  if (entry.HasMember(kSkew)) {
    prior->skew.is_set = true;
    prior->skew.value[0] = entry[kSkew].GetDouble();
  }

  // Get radial distortion coeffs.
  if (entry.HasMember(kRadialDistortionCoeffs) &&
      entry[kRadialDistortionCoeffs].IsArray()) {
    const int num_dist_coeffs = std::min(
        static_cast<int>(entry[kRadialDistortionCoeffs].Size()), 2);
    for (int i = 0; i < num_dist_coeffs; ++i) {
      prior->radial_distortion.value[i] =
          entry[kRadialDistortionCoeffs][i].GetDouble();
    }
    prior->radial_distortion.is_set =
        !entry[kRadialDistortionCoeffs].Empty();
  }

  return true;
}

bool ExtractCameraIntrinsicsPrior(const rapidjson::Value& entry,
                                  std::string* view_name,
                                  CameraIntrinsicsPrior* prior) {
  // Get the view name.
  if (!entry.HasMember(kImageName)) {
    LOG(ERROR) << "Could not find the image name.";
    return false;
  }
  *view_name = entry[kImageName].GetString();
  VLOG(3) << "Loading camera intrinsics prior for image name: " << *view_name;

  // Get the camera type.
  std::string camera_type_str;
  if (!entry.HasMember(kCameraType)) {
    LOG(WARNING) << "Unknown camera for view: " << *view_name
                 << ". Setting to PINHOLE.";
    camera_type_str = kPinholeType;
  } else {
    camera_type_str = entry[kCameraType].GetString();
  }
  // Check if the camera type exists.
  const CameraIntrinsicsModelType camera_type =
      StringToCameraIntrinsicsModelType(camera_type_str);

  switch (camera_type) {
    case CameraIntrinsicsModelType::PINHOLE:
      return ExtractPinholeCamera(entry, prior);
    case CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL:
    case CameraIntrinsicsModelType::FISHEYE:
    case CameraIntrinsicsModelType::FOV:
    case CameraIntrinsicsModelType::DIVISION_UNDISTORTION:
      LOG(FATAL) << "Camera type not supported yet.";
    default:
      LOG(ERROR) << "Invalid camera type.";
      return false;
  };

  return true;
}

}  // namespace

bool ExtractCameraIntrinsicPriorsFromJson(
    const char* json_str,
    std::unordered_map<std::string, CameraIntrinsicsPrior>* view_to_priors) {
  using rapidjson::Document;
  using rapidjson::SizeType;
  using rapidjson::Value;

  Document json;
  json.Parse(json_str);

  if (!json.HasMember(kPriorsEntry) || !json[kPriorsEntry].IsArray()) {
    LOG(ERROR) << "Expected \"priors\" array entry in JSON.";
    return false;
  }

  const Value& entries = json[kPriorsEntry];
  std::string view_name;
  for (SizeType i = 0; i < entries.Size(); ++i) {
    CameraIntrinsicsPrior prior;
    if (!entries[i].HasMember(kCameraIntrinsicsPrior) ||
        !ExtractCameraIntrinsicsPrior(entries[i][kCameraIntrinsicsPrior],
                                      &view_name,
                                      &prior)) {
      LOG(WARNING) << "Could not parse entry at position: " << i;
      continue;
    }
    // Add to the map.
    (*view_to_priors)[view_name] = prior;
  }

  return true;
}

bool ReadCalibration(const std::string& calibration_file,
                     std::unordered_map<std::string, CameraIntrinsicsPrior>*
                         camera_intrinsics_priors) {
  FILE* file = fopen(calibration_file.c_str(), "rb");
  if (file == nullptr) {
    LOG(ERROR) << "Cannot read the list file from " << calibration_file;
    return false;
  }

  // Get the size of the file.
  fseek(file, 0, SEEK_END);
  const size_t buffer_size = ftell(file);
  fseek(file, 0, SEEK_SET);  

  // Allocate a buffer
  char* file_buffer = new char[buffer_size];

  // Read the whole file.
  CHECK_LE(fread(file_buffer, buffer_size, 1, file), buffer_size);
  fclose(file);

  const bool json_parsed =
      ExtractCameraIntrinsicPriorsFromJson(file_buffer,
                                           camera_intrinsics_priors);

  delete [] file_buffer;
  file_buffer = nullptr;

  return json_parsed;
  
  // std::ifstream ifs(calibration_file.c_str(), std::ios::in);
  // if (!ifs.is_open()) {
  //   LOG(ERROR) << "Cannot read the list file from " << calibration_file;
  //   return false;
  // }

  // while (!ifs.eof()) {
  //   // Read in the filename.
  //   std::string filename;
  //   ifs >> filename;
  //   if (filename.length() == 0) {
  //     break;
  //   }

  //   // Read camera_intrinsics_prior.
  //   CameraIntrinsicsPrior temp_camera_intrinsics_prior;
  //   temp_camera_intrinsics_prior.focal_length.is_set = true;
  //   ifs >> temp_camera_intrinsics_prior.focal_length.value[0];

  //   temp_camera_intrinsics_prior.principal_point.is_set = true;
  //   ifs >> temp_camera_intrinsics_prior.principal_point.value[0];
  //   ifs >> temp_camera_intrinsics_prior.principal_point.value[1];
  //   temp_camera_intrinsics_prior.image_width =
  //       2.0 * temp_camera_intrinsics_prior.principal_point.value[0];
  //   temp_camera_intrinsics_prior.image_height =
  //       2.0 * temp_camera_intrinsics_prior.principal_point.value[1];

  //   temp_camera_intrinsics_prior.aspect_ratio.is_set = true;
  //   ifs >> temp_camera_intrinsics_prior.aspect_ratio.value[0];

  //   temp_camera_intrinsics_prior.skew.is_set = true;
  //   ifs >> temp_camera_intrinsics_prior.skew.value[0];

  //   temp_camera_intrinsics_prior.radial_distortion.is_set = true;
  //   ifs >> temp_camera_intrinsics_prior.radial_distortion.value[0];
  //   ifs >> temp_camera_intrinsics_prior.radial_distortion.value[1];

  //   (*camera_intrinsics_prior)[filename] = temp_camera_intrinsics_prior;
  // }

  return true;
}

}  // namespace theia
