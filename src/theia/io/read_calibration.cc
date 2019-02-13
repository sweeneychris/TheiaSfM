// Copyright (C) 2019 The Regents of the University of California (Regents).
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
#include <vector>

#include <rapidjson/document.h>
#include <stlplus3/file_system.hpp>

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
static const char* kPrincipalPoint = "principal_point";
static const char* kAspectRatio = "aspect_ratio";
static const char* kSkew = "skew";
static const char* kRadialDistortionCoeffs = "radial_distortion_coeffs";
static const char* kTangentialDistortionCoeffs = "tangential_distortion_coeffs";
static const char* kPosition = "position";
static const char* kOrientation = "orientation";
static const char* kLatitude = "latitude";
static const char* kLongitude = "longitude";
static const char* kAltitude = "altitude";

bool ExtractPriorParameters(const rapidjson::Value& entry,
                            CameraIntrinsicsPrior* prior) {
  // Get the focal length.
  if (entry.HasMember(kFocalLength)) {
    prior->focal_length.is_set = true;
    prior->focal_length.value[0] = entry[kFocalLength].GetDouble();
  }

  // Get the principal points.
  if (entry.HasMember(kPrincipalPoint) && entry[kPrincipalPoint].IsArray()) {
    const int num_entries =
        std::min(static_cast<int>(entry[kPrincipalPoint].Size()), 2);

    bool all_doubles = true;
    for (int i = 0; i < num_entries; ++i) {
      bool is_double =
          entry[kPrincipalPoint][i].IsDouble() ||
          entry[kPrincipalPoint][i].IsInt();
      if (is_double) {
        prior->principal_point.value[i] = entry[kPrincipalPoint][i].GetDouble();
      }
      all_doubles = all_doubles && is_double;
    }
    prior->principal_point.is_set =
        !entry[kPrincipalPoint].Empty() && all_doubles;

    if (prior->principal_point.is_set) {
      prior->image_width =
          static_cast<int>(2 * prior->principal_point.value[0]);
      prior->image_height =
          static_cast<int>(2 * prior->principal_point.value[1]);
    }
  }

  // Get width.
  if (entry.HasMember(kImageWidth) && entry[kImageWidth].IsInt()) {
    prior->image_width = entry[kImageWidth].GetInt();
  }

  // Get height.
  if (entry.HasMember(kImageHeight) && entry[kImageHeight].IsInt()) {
    prior->image_height = entry[kImageHeight].GetInt();
  }

  // Get aspect ratio.
  if (entry.HasMember(kAspectRatio) &&
      (entry[kAspectRatio].IsDouble() || entry[kAspectRatio].IsInt())) {
    prior->aspect_ratio.is_set = true;
    prior->aspect_ratio.value[0] = entry[kAspectRatio].GetDouble();
  }

  // Get skew.
  if (entry.HasMember(kSkew) &&
      (entry[kSkew].IsDouble() || entry[kSkew].IsInt())) {
    prior->skew.is_set = true;
    prior->skew.value[0] = entry[kSkew].GetDouble();
  }

  // Get radial distortion coeffs.
  if (entry.HasMember(kRadialDistortionCoeffs) &&
      entry[kRadialDistortionCoeffs].IsArray()) {
    const int num_dist_coeffs = std::min(
        static_cast<int>(entry[kRadialDistortionCoeffs].Size()), 4);
    bool all_doubles = true;
    for (int i = 0; i < num_dist_coeffs; ++i) {
      bool is_double =
          entry[kRadialDistortionCoeffs][i].IsDouble() ||
          entry[kRadialDistortionCoeffs][i].IsInt();
      if (is_double) {
        prior->radial_distortion.value[i] =
            entry[kRadialDistortionCoeffs][i].GetDouble();
      }
      all_doubles = all_doubles && is_double;
    }
    prior->radial_distortion.is_set =
        !entry[kRadialDistortionCoeffs].Empty() && all_doubles;
  }

  // Get tangential distortion coeffs.
  if (entry.HasMember(kTangentialDistortionCoeffs) &&
      entry[kTangentialDistortionCoeffs].IsArray()) {
    const int num_dist_coeffs = std::min(
        static_cast<int>(entry[kTangentialDistortionCoeffs].Size()), 2);
    bool all_doubles = true;
    for (int i = 0; i < num_dist_coeffs; ++i) {
      bool is_double =
          entry[kTangentialDistortionCoeffs][i].IsDouble() ||
          entry[kTangentialDistortionCoeffs][i].IsInt();
      if (is_double) {
        prior->tangential_distortion.value[i] =
            entry[kTangentialDistortionCoeffs][i].GetDouble();
      }
      all_doubles = all_doubles && is_double;
    }
    prior->tangential_distortion.is_set =
        !entry[kTangentialDistortionCoeffs].Empty() && all_doubles;
  }

  // Get position.
  if (entry.HasMember(kPosition) && entry[kPosition].IsArray()) {
    const int num_entries = std::min(
        static_cast<int>(entry[kPosition].Size()), 3);
    bool all_doubles = true;
    for (int i = 0; i < num_entries; ++i) {
      bool is_double =
          entry[kPosition][i].IsDouble() ||
          entry[kPosition][i].IsInt();
      if (is_double) {
        prior->position.value[i] = entry[kPosition][i].GetDouble();
      }
      all_doubles = all_doubles && is_double;
    }
    prior->position.is_set = !entry[kPosition].Empty() && all_doubles;
  }

  // Get orientation using Angle-Axis.
  if (entry.HasMember(kOrientation) && entry[kOrientation].IsArray()) {
    const int num_entries = std::min(
        static_cast<int>(entry[kOrientation].Size()), 3);
    bool all_doubles = true;
    for (int i = 0; i < num_entries; ++i) {
      bool is_double =
          entry[kOrientation][i].IsDouble() || entry[kOrientation][i].IsInt();
      if (is_double) {
        prior->orientation.value[i] = entry[kOrientation][i].GetDouble();
      }
      all_doubles = all_doubles && is_double;
    }
    prior->orientation.is_set = !entry[kOrientation].Empty() && all_doubles;
  }

  // Get GPS priors.
  if (entry.HasMember(kLatitude) &&
      (entry[kLatitude].IsDouble() || entry[kLatitude].IsInt())) {
    prior->latitude.value[0] = entry[kLatitude].GetDouble();
    prior->latitude.is_set = true;
  }

  if (entry.HasMember(kLongitude) &&
      (entry[kLongitude].IsDouble() || entry[kLongitude].IsInt())) {
    prior->longitude.value[0] = entry[kLongitude].GetDouble();
    prior->longitude.is_set = true;
  }

  
  if (entry.HasMember(kAltitude) &&
      (entry[kAltitude].IsDouble() || entry[kAltitude].IsInt())) {
    prior->altitude.value[0] = entry[kAltitude].GetDouble();
    prior->altitude.is_set = true;
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
    VLOG(3) << "Camera type: [" << camera_type_str << "]";
  }

  // Get the camera type. This will verify if that the camera type is valid.
  StringToCameraIntrinsicsModelType(camera_type_str);
  ExtractPriorParameters(entry, prior);

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
  // Get the size of the file.
  const size_t buffer_size = stlplus::file_size(calibration_file);

  // Allocate a buffer
  std::vector<char> file_buffer(buffer_size, 0);

  // Open and read the whole file.
  FILE* file = fopen(calibration_file.c_str(), "rb");
  if (file == nullptr) {
    LOG(ERROR) << "Cannot read file: " << calibration_file;
    return false;
  }

  // Read the whole file.
  CHECK_EQ(fread(file_buffer.data(), sizeof(file_buffer[0]), buffer_size, file),
           buffer_size);
  fclose(file);

  // Remove spurious chars after the closing curly brace.
  std::string file_content(file_buffer.begin(), file_buffer.end());
  const size_t last_curly_idx = file_content.rfind('}');
  if (last_curly_idx == std::string::npos) {
    LOG(ERROR) << "Could not fid a proper JSON file: " << calibration_file;
    return false;
  }

  file_content.resize(last_curly_idx + 1);

  const bool json_parsed =
      ExtractCameraIntrinsicPriorsFromJson(file_content.c_str(),
                                           camera_intrinsics_priors);

  return json_parsed;
}

}  // namespace theia
