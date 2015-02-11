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

#include "theia/io/write_matches.h"

#include <glog/logging.h>
#include <cstdlib>
#include <fstream>   // NOLINT
#include <iostream>  // NOLINT
#include <string>

#include "theia/matching/feature_matcher.h"
#include "theia/sfm/match_and_verify_features.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {
namespace {

void WriteCameraIntrinsicsPrior(const Prior& prior,
                                std::ofstream* matches_reader) {
  matches_reader->write(reinterpret_cast<const char*>(&prior.is_set),
                        sizeof(prior.is_set));
  matches_reader->write(reinterpret_cast<const char*>(&prior.value),
                        sizeof(prior.value));
}

// Write all image names.
void WriteView(const std::string& view_name,
               const CameraIntrinsicsPrior& camera_intrinsics_prior,
               std::ofstream* matches_writer) {
  const uint32_t name_length = view_name.length();
  matches_writer->write(reinterpret_cast<const char*>(&name_length),
                        sizeof(name_length));
  matches_writer->write(reinterpret_cast<const char*>(view_name.c_str()),
                        name_length);

  // Write image size.
  matches_writer->write(
      reinterpret_cast<const char*>(&camera_intrinsics_prior.image_width),
      sizeof(camera_intrinsics_prior.image_width));
  matches_writer->write(
      reinterpret_cast<const char*>(&camera_intrinsics_prior.image_height),
      sizeof(camera_intrinsics_prior.image_height));

  // Write view camera_intrinsics_prior.
  WriteCameraIntrinsicsPrior(camera_intrinsics_prior.focal_length,
                             matches_writer);
  WriteCameraIntrinsicsPrior(camera_intrinsics_prior.principal_point[0],
                             matches_writer);
  WriteCameraIntrinsicsPrior(camera_intrinsics_prior.principal_point[1],
                             matches_writer);
}

// Writes the two image indices.
void WriteImagePairIndices(const int image1_index, const int image2_index,
                           std::ofstream* matches_writer) {
  uint32_t image_index = image1_index;
  matches_writer->write(reinterpret_cast<char*>(&image_index),
                        sizeof(image_index));
  image_index = image2_index;
  matches_writer->write(reinterpret_cast<char*>(&image_index),
                        sizeof(image_index));
}

// Write a two view info struct from the matches files.
void WriteTwoViewInfo(const TwoViewInfo& info, std::ofstream* matches_writer) {
  // Write focal lengths.
  matches_writer->write(reinterpret_cast<const char*>(&info.focal_length_1),
                        sizeof(info.focal_length_1));
  matches_writer->write(reinterpret_cast<const char*>(&info.focal_length_2),
                        sizeof(info.focal_length_2));
  // Write relative position and rotation.
  matches_writer->write(reinterpret_cast<const char*>(info.position_2.data()),
                        sizeof(info.position_2));
  matches_writer->write(reinterpret_cast<const char*>(info.rotation_2.data()),
                        sizeof(info.rotation_2));

  // Write number of geometrically verified features.
  matches_writer->write(reinterpret_cast<const char*>(&info.num_verified_matches),
                        sizeof(info.num_verified_matches));
}

void WriteFeatureMatches(const std::vector<FeatureCorrespondence>& matches,
                         std::ofstream* matches_writer) {
  const uint32_t num_feature_matches = matches.size();
  matches_writer->write(reinterpret_cast<const char*>(&num_feature_matches),
                        sizeof(num_feature_matches));

  for (int i = 0; i < num_feature_matches; i++) {
    matches_writer->write(
        reinterpret_cast<const char*>(matches[i].feature1.data()),
        sizeof(matches[i].feature1));
    matches_writer->write(
        reinterpret_cast<const char*>(matches[i].feature2.data()),
        sizeof(matches[i].feature2));
  }
}

}  // namespace

bool WriteMatchesAndGeometry(
    const std::string& matches_file,
    const std::vector<std::string>& view_names,
    const std::vector<CameraIntrinsicsPrior>& camera_intrinsics_prior,
    const std::vector<ImagePairMatch>& matches) {
  CHECK_EQ(view_names.size(), camera_intrinsics_prior.size());

  // Return false if the file cannot be opened for writing.
  std::ofstream matches_writer(matches_file, std::ios::out | std::ios::binary);
  if (!matches_writer.is_open()) {
    LOG(ERROR) << "Could not open the matches file: " << matches_file
               << " for writing.";
    return false;
  }

  // Write all view data.
  const uint32_t num_views = view_names.size();
  matches_writer.write(reinterpret_cast<const char*>(&num_views),
                       sizeof(num_views));
  for (int i = 0; i < num_views; i++) {
    WriteView(view_names[i], camera_intrinsics_prior[i], &matches_writer);
  }

  // Write number of image pair matches.
  const uint64_t num_image_matches = matches.size();
  matches_writer.write(reinterpret_cast<const char*>(&num_image_matches),
                       sizeof(num_image_matches));

  // Write all image matches.
  for (const ImagePairMatch& match : matches) {
    WriteImagePairIndices(match.image1_index,
                          match.image2_index,
                          &matches_writer);
    WriteTwoViewInfo(match.twoview_info, &matches_writer);
    WriteFeatureMatches(match.correspondences, &matches_writer);
  }

  return true;
}

}  // namespace theia
