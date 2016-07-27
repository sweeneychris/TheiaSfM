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

#include <easyexif/exif.h>
#include <glog/logging.h>

#include <algorithm>
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
  static const double kMinFocalLength = 1e-2;

  // Read in the EXIF information.
  std::streampos file_size;
  std::ifstream file(image_file, std::ios::binary);
  if (!file.is_open()) {
    return false;
  }

  // Get the file size size:
  file.seekg(0, std::ios::end);
  file_size = file.tellg();
  file.seekg(0, std::ios::beg);

  // Read the data:
  std::vector<unsigned char> jpeg_data(file_size);
  file.read(reinterpret_cast<char*>(&jpeg_data[0]), file_size);

  // TODO(cmsweeney): The exif parser returns a success code. We may want to log
  // the success code.
  easyexif::EXIFInfo exif_parser;
  exif_parser.parseFrom(jpeg_data.data(), file_size);
  file.close();

  // Set the image dimensions.
  const FloatImage image(image_file);
  camera_intrinsics_prior->image_width = image.Width();
  camera_intrinsics_prior->image_height = image.Height();

  // Set principal point.
  camera_intrinsics_prior->principal_point.is_set = true;
  camera_intrinsics_prior->principal_point.value[0] =
      camera_intrinsics_prior->image_width / 2.0;
  camera_intrinsics_prior->principal_point.value[1] =
      camera_intrinsics_prior->image_height / 2.0;

  // If the exif focal length (in mm) is not set, then the focal length in
  // pixels cannot be set.
  if (exif_parser.FocalLength < kMinFocalLength) {
    return true;
  }

  // Attempt to set the focal length from the plane resolution, then try the
  // sensor width database if that fails.
  if (exif_parser.ImageWidth > 0 &&
      exif_parser.LensInfo.FocalPlaneXResolution > 0 &&
      exif_parser.LensInfo.FocalPlaneResolutionUnit > 1 &&
      exif_parser.LensInfo.FocalPlaneResolutionUnit <= 5) {
    SetFocalLengthFromExif(exif_parser,
                           image.Width(),
                           image.Height(),
                           camera_intrinsics_prior);
  } else {
    SetFocalLengthFromSensorDatabase(
        exif_parser,
        std::max(camera_intrinsics_prior->image_width,
                 camera_intrinsics_prior->image_height),
        camera_intrinsics_prior);
  }

  if (camera_intrinsics_prior->focal_length.value[0] > kMinFocalLength) {
    camera_intrinsics_prior->focal_length.is_set = true;
  }

  // Set GPS latitude, longitude, and altitude.
  if (exif_parser.GeoLocation.Latitude != 0) {
    camera_intrinsics_prior->latitude.is_set = true;
    camera_intrinsics_prior->latitude.value[0] =
        exif_parser.GeoLocation.Latitude;
  }
  if (exif_parser.GeoLocation.Longitude != 0) {
    camera_intrinsics_prior->longitude.is_set = true;
    camera_intrinsics_prior->longitude.value[0] =
        exif_parser.GeoLocation.Longitude;
  }
  if (exif_parser.GeoLocation.Altitude != 0) {
    camera_intrinsics_prior->altitude.is_set = true;
    camera_intrinsics_prior->altitude.value[0] =
        exif_parser.GeoLocation.Altitude;
  }

  return true;
}

void ExifReader::SetFocalLengthFromExif(
    const easyexif::EXIFInfo& exif_parser,
    const double image_width,
    const double image_height,
    CameraIntrinsicsPrior* camera_intrinsics_prior) const {
  // CCD resolution is the pixels per unit resolution of the CCD.
  double ccd_resolution_units = 1.0;
  switch (exif_parser.LensInfo.FocalPlaneResolutionUnit) {
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
      break;
  }

  const double ccd_width =
      exif_parser.ImageWidth /
      (exif_parser.LensInfo.FocalPlaneXResolution / ccd_resolution_units);
  const double ccd_height =
      exif_parser.ImageHeight /
      (exif_parser.LensInfo.FocalPlaneYResolution / ccd_resolution_units);

  const double focal_length_x =
      exif_parser.FocalLength * image_width / ccd_width;
  const double focal_length_y =
      exif_parser.FocalLength * image_height / ccd_height;

  // Normalize for the image size in case the original size is different
  // than the current size.
  camera_intrinsics_prior->focal_length.value[0] =
      (focal_length_x + focal_length_y) / 2.0;
  camera_intrinsics_prior->focal_length.is_set = true;
}

void ExifReader::SetFocalLengthFromSensorDatabase(
    const easyexif::EXIFInfo& exif_parser,
    const double max_image_dimension,
    CameraIntrinsicsPrior* camera_intrinsics_prior) const {
  const std::string make = ToLowercase(exif_parser.Make);
  const std::string model = ToLowercase(exif_parser.Model);

  // First, try to look up just the model.
  const std::string make_model = make + " " + model;
  double sensor_width = 0;
  if (ContainsKey(sensor_width_database_, model)) {
    sensor_width = FindOrDie(sensor_width_database_, model);
  } else if (ContainsKey(sensor_width_database_, make_model)) {
    sensor_width = FindOrDie(sensor_width_database_, make_model);
  }

  if (sensor_width != 0) {
    camera_intrinsics_prior->focal_length.value[0] =
        max_image_dimension * exif_parser.FocalLength / sensor_width;
  }
}

}  // namespace theia
