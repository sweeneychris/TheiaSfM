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

#include "theia/sfm/exif_reader.h"

#include <glog/logging.h>

#include <algorithm>
#include <cmath>
#include <fstream>  // NOLINT
#include <iostream>  // NOLINT
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "theia/image/image.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/util/map_util.h"

namespace theia {
namespace {

void RemoveLeadingTrailingSpaces(std::string* str) {
  size_t p = str->find_first_not_of(" \t");
  str->erase(0, p);

  p = str->find_last_not_of(" \t");
  if (std::string::npos != p)
    str->erase(p + 1);
}

std::string ToLowercase(const std::string& str) {
  std::string str2 = str;
  RemoveLeadingTrailingSpaces(&str2);
  std::transform(str2.begin(), str2.end(), str2.begin(), ::tolower);
  if (!str2.empty() &&
      (str2[str2.size() - 1] == '\r' || str2[str2.size() - 1] == '\0')) {
    str2.erase(str2.size() - 1);
  }
  return str2;
}

}  // namespace

std::vector<std::string> SplitString(const std::string &s, const char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

ExifReader::ExifReader() {
  LoadSensorWidthDatabase();
}

void ExifReader::LoadSensorWidthDatabase() {
  const std::string sensor_width_file =
      std::string(THEIA_DATA_DIR) + "/camera_sensor_database.txt";

  std::ifstream ifs(sensor_width_file.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(FATAL) << "Cannot read the sensor width file from "
               << sensor_width_file;
  }

  while (!ifs.eof()) {
    // Read in the filename.
    std::string line;
    std::getline(ifs, line);
    if (line.size() == 0) {
      break;
    }

    const auto& tokens = SplitString(line, ';');
    CHECK_EQ(tokens.size(), 3);

    const std::string make = ToLowercase(tokens[0]);
    const std::string model = ToLowercase(tokens[1]);
    const double camera_sensor_width = stod(tokens[2]);

    // In the database, the model includes the make.
    InsertOrDie(&sensor_width_database_, model, camera_sensor_width);
  }
}

bool ExifReader::ExtractEXIFMetadata(
    const std::string& image_file,
    CameraIntrinsicsPrior* camera_intrinsics_prior) const {
  CHECK_NOTNULL(camera_intrinsics_prior);

  OpenImageIO::ImageBuf image(image_file);
  OpenImageIO::ImageSpec image_spec = image.spec();

  // Set the image dimensions.
  camera_intrinsics_prior->image_width = image_spec.width;
  camera_intrinsics_prior->image_height = image_spec.height;

  // Set principal point.
  camera_intrinsics_prior->principal_point.is_set = true;
  camera_intrinsics_prior->principal_point.value[0] =
      camera_intrinsics_prior->image_width / 2.0;
  camera_intrinsics_prior->principal_point.value[1] =
      camera_intrinsics_prior->image_height / 2.0;

  // Attempt to set the focal length from the plane resolution, then try the
  // sensor width database if that fails.
  if (!SetFocalLengthFromExif(image_spec, camera_intrinsics_prior) &&
      !SetFocalLengthFromSensorDatabase(image_spec, camera_intrinsics_prior)) {
      return true;
  }

  // Make sure the focal length value is sane. If so, indicate the prior as
  // being set.
  const double focal_length = camera_intrinsics_prior->focal_length.value[0];
  camera_intrinsics_prior->focal_length.is_set = true;

  // Set GPS latitude.
  const OpenImageIO::ImageIOParameter* latitude =
      image_spec.find_attribute("GPS:Latitude");
  if (latitude != nullptr) {
    camera_intrinsics_prior->latitude.is_set = true;
    const float* latitude_val =
        reinterpret_cast<const float*>(latitude->data());
    camera_intrinsics_prior->latitude.value[0] =
        latitude_val[0] + latitude_val[1] / 60.0 + latitude_val[2] / 3600.0;

    // Adjust the sign of the latitude depending on if the coordinates given
    // were north or south.
    const std::string north_or_south =
        image_spec.get_string_attribute("GPS:LatitudeRef");
    if (north_or_south == "S") {
      camera_intrinsics_prior->longitude.value[0] *= -1.0;
    }
  }

  // Set GPS longitude.
  const OpenImageIO::ImageIOParameter* longitude =
      image_spec.find_attribute("GPS:Longitude");
  if (longitude != nullptr) {
    camera_intrinsics_prior->longitude.is_set = true;
    const float* longitude_val =
        reinterpret_cast<const float*>(longitude->data());
    camera_intrinsics_prior->longitude.value[0] =
        longitude_val[0] + longitude_val[1] / 60.0 + longitude_val[2] / 3600.0;

    // Adjust the sign of the longitude depending on if the coordinates given
    // were east or west.
    const std::string east_or_west =
        image_spec.get_string_attribute("GPS:LongitudeRef");
    if (east_or_west == "W") {
      camera_intrinsics_prior->longitude.value[0] *= -1.0;
    }
  }


  // Set GSP altitude.
  const OpenImageIO::ImageIOParameter* altitude =
      image_spec.find_attribute("GPS:Altitude");
  if (altitude != nullptr) {
    camera_intrinsics_prior->altitude.is_set = true;
    camera_intrinsics_prior->altitude.value[0] =
        image_spec.get_float_attribute("GPS:Altitude");
  }

  return true;
}

bool ExifReader::SetFocalLengthFromExif(
    const OpenImageIO::ImageSpec& image_spec,
    CameraIntrinsicsPrior* camera_intrinsics_prior) const {
  static const float kMinFocalLength = 1e-2;

  const float exif_focal_length =
    image_spec.get_float_attribute("Exif:FocalLength", kMinFocalLength);
  const float focal_plane_x_resolution =
      image_spec.get_float_attribute("Exif:FocalPlaneXResolution");
  const float focal_plane_y_resolution =
      image_spec.get_float_attribute("Exif:FocalPlaneYResolution");
  const int focal_plane_resolution_unit =
      image_spec.get_int_attribute("Exif:FocalPlaneResolutionUnit");

  // Make sure the values are sane.
  if (exif_focal_length <= kMinFocalLength || focal_plane_x_resolution <= 0.0 ||
      focal_plane_y_resolution <= 0.0) {
    return false;
  }

  // CCD resolution is the pixels per unit resolution of the CCD.
  double ccd_resolution_units = 1.0;
  switch (focal_plane_resolution_unit) {
    case 2:
      // Convert inches to mm.
      ccd_resolution_units = 25.4;
      break;
    case 3:
      // Convert centimeters to mm.
      ccd_resolution_units = 10.0;
      break;
    case 4:
      // Already in mm.
      break;
    case 5:
      // Convert micrometers to mm.
      ccd_resolution_units = 1.0 / 1000.0;
      break;
    default:
      return false;
      break;
  }

  // Get the ccd dimensions in mm.
  const int exif_width = image_spec.get_int_attribute("Exif:PixelXDimension");
  const int exif_height = image_spec.get_int_attribute("Exif:PixelYDimension");
  const double ccd_width =
      exif_width / (focal_plane_x_resolution / ccd_resolution_units);
  const double ccd_height =
      exif_height / (focal_plane_y_resolution / ccd_resolution_units);

  const double focal_length_x =
      exif_focal_length * image_spec.width / ccd_width;
  const double focal_length_y =
      exif_focal_length * image_spec.height / ccd_height;

  // Normalize for the image size in case the original size is different
  // than the current size.
  const double focal_length = (focal_length_x + focal_length_y) / 2.0;
  camera_intrinsics_prior->focal_length.value[0] = focal_length;
  return std::isfinite(focal_length) && focal_length > 0.0;
}

bool ExifReader::SetFocalLengthFromSensorDatabase(
    const OpenImageIO::ImageSpec& image_spec,
    CameraIntrinsicsPrior* camera_intrinsics_prior) const {
  const int max_image_dimension = std::max(image_spec.width, image_spec.height);
  const float exif_focal_length =
      image_spec.get_float_attribute("Exif:FocalLength");

  const std::string camera_make = image_spec.get_string_attribute("Make");
  const std::string camera_model = image_spec.get_string_attribute("Model");
  // First, try to look up just the model.
  const std::string make_model =
      ToLowercase(camera_make) + " " + ToLowercase(camera_model);
  double sensor_width = 0;
  if (ContainsKey(sensor_width_database_, camera_model)) {
    sensor_width = FindOrDie(sensor_width_database_, camera_model);
  } else if (ContainsKey(sensor_width_database_, make_model)) {
    sensor_width = FindOrDie(sensor_width_database_, make_model);
  }

  if (sensor_width == 0) {
    return false;
  }

  const double focal_length =
      max_image_dimension * exif_focal_length / sensor_width;
  camera_intrinsics_prior->focal_length.value[0] = focal_length;
  return std::isfinite(focal_length) && focal_length > 0.0;
}

}  // namespace theia
