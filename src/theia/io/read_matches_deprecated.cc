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

#include "theia/io/read_matches_deprecated.h"

#include <glog/logging.h>
#include <stdint.h>
#include <cstdlib>
#include <fstream>   // NOLINT
#include <iostream>  // NOLINT
#include <string>

#include "theia/matching/feature_matcher.h"
#include "theia/sfm/match_and_verify_features.h"
#include "theia/sfm/twoview_info.h"

namespace theia {
namespace {

// Reads a camera_intrinsics_prior field.
void ReadCameraIntrinsicsPrior(std::ifstream* matches_reader,
                               Prior* prior) {
  matches_reader->read(reinterpret_cast<char*>(&prior->is_set),
                       sizeof(prior->is_set));
  matches_reader->read(reinterpret_cast<char*>(&prior->value),
                       sizeof(prior->value));
}

// Read the data for a particular view.
void ReadView(std::ifstream* matches_reader,
              std::string* view_name,
              CameraIntrinsicsPrior* camera_intrinsics_prior) {
  // Read the view name.
  uint32_t name_length;
  matches_reader->read(reinterpret_cast<char*>(&name_length),
                       sizeof(name_length));
  char* buffer = new char[name_length + 1];
  matches_reader->read(buffer, name_length);
  buffer[name_length] = '\0';
  *view_name = std::string(buffer, name_length);
  delete[] buffer;

  // Read image size.
  matches_reader->read(
      reinterpret_cast<char*>(&camera_intrinsics_prior->image_width),
      sizeof(camera_intrinsics_prior->image_width));
  matches_reader->read(
      reinterpret_cast<char*>(&camera_intrinsics_prior->image_height),
      sizeof(camera_intrinsics_prior->image_height));

  // Read view camera_intrinsics_prior.
  ReadCameraIntrinsicsPrior(matches_reader,
                            &camera_intrinsics_prior->focal_length);
  ReadCameraIntrinsicsPrior(matches_reader,
                            &camera_intrinsics_prior->principal_point[0]);
  ReadCameraIntrinsicsPrior(matches_reader,
                            &camera_intrinsics_prior->principal_point[1]);
  ReadCameraIntrinsicsPrior(matches_reader,
                            &camera_intrinsics_prior->aspect_ratio);
  ReadCameraIntrinsicsPrior(matches_reader,
                            &camera_intrinsics_prior->skew);
  ReadCameraIntrinsicsPrior(matches_reader,
                            &camera_intrinsics_prior->radial_distortion[0]);
  ReadCameraIntrinsicsPrior(matches_reader,
                            &camera_intrinsics_prior->radial_distortion[1]);
}

// Reads the two image indices.
void ReadImagePairIndices(std::ifstream* matches_reader,
                          int* image1_index, int* image2_index) {
  uint32_t image_index;
  matches_reader->read(reinterpret_cast<char*>(&image_index),
                       sizeof(image_index));
  *image1_index = image_index;
  matches_reader->read(reinterpret_cast<char*>(&image_index),
                       sizeof(image_index));
  *image2_index = image_index;
}

// Read a two view info struct from the matches files.
void ReadTwoViewInfo(std::ifstream* matches_reader, TwoViewInfo* info) {
  // Read focal lengths.
  matches_reader->read(reinterpret_cast<char*>(&info->focal_length_1),
                       sizeof(info->focal_length_1));
  matches_reader->read(reinterpret_cast<char*>(&info->focal_length_2),
                       sizeof(info->focal_length_2));
  // Read relative position and rotation.
  matches_reader->read(reinterpret_cast<char*>(info->position_2.data()),
                       sizeof(info->position_2));
  matches_reader->read(reinterpret_cast<char*>(info->rotation_2.data()),
                       sizeof(info->rotation_2));

  // Read number of geometrically verified features.
  matches_reader->read(reinterpret_cast<char*>(&info->num_verified_matches),
                       sizeof(info->num_verified_matches));
}

void ReadFeatureMatches(std::ifstream* matches_reader,
                        std::vector<FeatureCorrespondence>* matches) {
  uint32_t num_feature_matches;
  matches_reader->read(reinterpret_cast<char*>(&num_feature_matches),
                       sizeof(num_feature_matches));
  matches->resize(num_feature_matches);
  for (int i = 0; i < num_feature_matches; i++) {
    matches_reader->read(reinterpret_cast<char*>((*matches)[i].feature1.data()),
                         sizeof((*matches)[i].feature1));
    matches_reader->read(reinterpret_cast<char*>((*matches)[i].feature2.data()),
                         sizeof((*matches)[i].feature2));
  }
}

}  // namespace

bool ReadMatchesAndGeometryDeprecated(
    const std::string& matches_file,
    std::vector<std::string>* view_names,
    std::vector<CameraIntrinsicsPrior>* camera_intrinsics_prior,
    std::vector<ImagePairMatch>* matches) {
  LOG(WARNING) << "You are using a deprecated version of the matches "
                  "reader. Proceed with caution.";

  CHECK_NOTNULL(view_names)->clear();
  CHECK_NOTNULL(camera_intrinsics_prior)->clear();
  CHECK_NOTNULL(matches)->clear();

  // Return false if the file cannot be opened.
  std::ifstream matches_reader(matches_file, std::ios::in | std::ios::binary);
  if (!matches_reader.is_open()) {
    LOG(ERROR) << "Could not open the matches file: " << matches_file
               << " for reading.";
    return false;
  }

  // Read all view information.
  uint32_t num_views;
  matches_reader.read(reinterpret_cast<char*>(&num_views), sizeof(num_views));
  view_names->resize(num_views);
  camera_intrinsics_prior->resize(num_views);
  for (int i = 0; i < num_views; i++) {
    ReadView(&matches_reader, &(*view_names)[i],
             &(*camera_intrinsics_prior)[i]);
  }

  // Read image pair matches.
  uint64_t num_image_matches;
  matches_reader.read(reinterpret_cast<char*>(&num_image_matches),
                      sizeof(num_image_matches));
  matches->resize(num_image_matches);
  for (auto& match : *matches) {
    ReadImagePairIndices(&matches_reader,
                         &match.image1_index,
                         &match.image2_index);
    ReadTwoViewInfo(&matches_reader, &match.twoview_info);
    ReadFeatureMatches(&matches_reader, &match.correspondences);
  }

  return true;
}

}  // namespace theia
