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
// Author: Chris Sweeney (sweeneychris@gmail.com)

#include "theia/matching/in_memory_features_and_matches_database.h"

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <cstdlib>
#include <fstream>  // NOLINT
#include <glog/logging.h>
#include <iostream>  // NOLINT
#include <mutex>     // NOLINT
#include <string>

#include "theia/matching/image_pair_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/util/map_util.h"

namespace theia {

bool InMemoryFeaturesAndMatchesDatabase::ContainsCameraIntrinsicsPrior(
    const std::string& image_name) {
  return ContainsKey(intrinsics_priors_, image_name);
}

// Get/set the features for the image.
CameraIntrinsicsPrior
InMemoryFeaturesAndMatchesDatabase::GetCameraIntrinsicsPrior(
    const std::string& image_name) {
  return FindOrDie(intrinsics_priors_, image_name);
}

// Set the features for the image.
void InMemoryFeaturesAndMatchesDatabase::PutCameraIntrinsicsPrior(
    const std::string& image_name, const CameraIntrinsicsPrior& intrinsics) {
  intrinsics_priors_[image_name] = intrinsics;
}

// Supply an iterator to iterate over the priors.
std::vector<std::string>
InMemoryFeaturesAndMatchesDatabase::ImageNamesOfCameraIntrinsicsPriors() {
  std::vector<std::string> image_names;
  image_names.reserve(intrinsics_priors_.size());
  for (const auto& intrinsics : intrinsics_priors_) {
    image_names.push_back(intrinsics.first);
  }
  return image_names;
}

size_t InMemoryFeaturesAndMatchesDatabase::NumCameraIntrinsicsPrior() {
  return intrinsics_priors_.size();
}

bool InMemoryFeaturesAndMatchesDatabase::ContainsFeatures(
    const std::string& image_name) {
  return ContainsKey(features_, image_name);
}

// Get/set the features for the image.
KeypointsAndDescriptors InMemoryFeaturesAndMatchesDatabase::GetFeatures(
    const std::string& image_name) {
  return FindOrDie(features_, image_name);
}

// Set the features for the image.
void InMemoryFeaturesAndMatchesDatabase::PutFeatures(
    const std::string& image_name, const KeypointsAndDescriptors& features) {
  features_[image_name] = features;
}

std::vector<std::string>
InMemoryFeaturesAndMatchesDatabase::ImageNamesOfFeatures() {
  std::vector<std::string> features_keys;
  features_keys.reserve(features_.size());
  for (const auto& features : features_) {
    features_keys.push_back(features.first);
  }
  return features_keys;
}

size_t InMemoryFeaturesAndMatchesDatabase::NumImages() {
  return features_.size();
}

// Get the image pair match for the images.
ImagePairMatch InMemoryFeaturesAndMatchesDatabase::GetImagePairMatch(
    const std::string& image_name1, const std::string& image_name2) {
  return FindOrDieNoPrint(matches_, std::make_pair(image_name1, image_name2));
}

// Set the image pair match for the images.
void InMemoryFeaturesAndMatchesDatabase::PutImagePairMatch(
    const std::string& image_name1,
    const std::string& image_name2,
    const ImagePairMatch& matches) {
  std::lock_guard<std::mutex> lock(mutex_);
  matches_[std::make_pair(image_name1, image_name2)] = matches;
}

std::vector<std::pair<std::string, std::string>>
InMemoryFeaturesAndMatchesDatabase::ImageNamesOfMatches() {
  std::vector<std::pair<std::string, std::string>> match_keys;
  match_keys.reserve(matches_.size());
  for (const auto& match : matches_) {
    match_keys.push_back(match.first);
  }
  return match_keys;
}

size_t InMemoryFeaturesAndMatchesDatabase::NumMatches() {
  return matches_.size();
}

bool InMemoryFeaturesAndMatchesDatabase::ReadFromFile(
    const std::string& filepath) {
  // Return false if the file cannot be opened.
  std::ifstream matches_reader(filepath, std::ios::in | std::ios::binary);
  if (!matches_reader.is_open()) {
    LOG(ERROR) << "Could not open the matches file: " << filepath
               << " for reading.";
    return false;
  }

  // Make sure that Cereal is able to finish executing before returning.
  std::vector<ImagePairMatch> matches;
  std::vector<std::string> view_names;
  std::vector<CameraIntrinsicsPrior> camera_intrinsics_prior;
  {
    cereal::PortableBinaryInputArchive input_archive(matches_reader);
    input_archive(view_names, camera_intrinsics_prior, matches);
  }
  CHECK_EQ(view_names.size(), camera_intrinsics_prior.size());

  matches_.reserve(matches.size());
  for (const auto& match : matches) {
    matches_[std::make_pair(match.image1, match.image2)] = match;
  }

  intrinsics_priors_.reserve(camera_intrinsics_prior.size());
  for (int i = 0; i < view_names.size(); i++) {
    intrinsics_priors_[view_names[i]] = camera_intrinsics_prior[i];
  }

  return true;
}
bool InMemoryFeaturesAndMatchesDatabase::WriteToFile(
    const std::string& filepath) {
  // Return false if the file cannot be opened for writing.
  std::ofstream matches_writer(filepath, std::ios::out | std::ios::binary);
  if (!matches_writer.is_open()) {
    LOG(ERROR) << "Could not open the matches file: " << filepath
               << " for writing.";
    return false;
  }

  // Make sure that Cereal is able to finish executing before returning.
  std::vector<ImagePairMatch> matches;
  matches.reserve(matches_.size());
  for (const auto& match : matches_) {
    matches.push_back(match.second);
  }

  std::vector<std::string> view_names;
  view_names.reserve(intrinsics_priors_.size());
  std::vector<CameraIntrinsicsPrior> camera_intrinsics_prior;
  camera_intrinsics_prior.reserve(intrinsics_priors_.size());
  for (const auto& prior : intrinsics_priors_) {
    view_names.push_back(prior.first);
    camera_intrinsics_prior.push_back(prior.second);
  }
  {
    cereal::PortableBinaryOutputArchive output_archive(matches_writer);
    output_archive(view_names, camera_intrinsics_prior, matches);
  }

  return true;
}

void InMemoryFeaturesAndMatchesDatabase::RemoveAllMatches() {
  matches_.clear();
}

}  // namespace theia
