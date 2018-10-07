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

#include "theia/matching/local_features_and_matches_database.h"

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <cstdlib>
#include <fstream>  // NOLINT
#include <glog/logging.h>
#include <iostream>  // NOLINT
#include <mutex>
#include <string>

#include "theia/io/read_keypoints_and_descriptors.h"
#include "theia/io/write_keypoints_and_descriptors.h"
#include "theia/matching/image_pair_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/util/filesystem.h"
#include "theia/util/map_util.h"
#include "theia/util/string.h"

namespace theia {
namespace {
inline std::string FeatureFilenameFromImage(const std::string& output_dir,
                                            const std::string& image) {
  return output_dir + image + ".features";
}
}  // namespace

LocalFeaturesAndMatchesDatabase::LocalFeaturesAndMatchesDatabase(
    const std::string& directory, const int max_cache_entries)
    : directory_(directory), features_cache_(nullptr) {
  AppendTrailingSlashIfNeeded(&directory_);

  // Determine if the directory for writing out feature exists. If not, try to
  // create it.
  if (!DirectoryExists(directory_)) {
    CHECK(CreateNewDirectory(directory_))
        << "Could not create the directory for storing features during "
           "matching: "
        << directory_;
  } else {
    // Load existing feature files into the database.
    std::vector<std::string> feature_files;
    CHECK(GetFilepathsFromWildcard(FeatureFilenameFromImage(directory_, "*"),
                                   &feature_files));
    image_names_.insert(feature_files.begin(), feature_files.end());
  }

  // Initialize the cache.
  std::function<KeypointsAndDescriptors(const std::string&)> fetch_images =
      std::bind(&LocalFeaturesAndMatchesDatabase::FetchImages,
                this,
                std::placeholders::_1);
  features_cache_.reset(new LRUFeatureCache(fetch_images, max_cache_entries));
}

LocalFeaturesAndMatchesDatabase::~LocalFeaturesAndMatchesDatabase() {}

bool LocalFeaturesAndMatchesDatabase::ContainsFeatures(
    const std::string& image_name) const {
  return ContainsKey(image_names_, image_name);
}
// Get/set the features for the image.
KeypointsAndDescriptors LocalFeaturesAndMatchesDatabase::GetFeatures(
    const std::string& image_name) {
  return features_cache_->Fetch(image_name);
}

// Set the features for the image.
void LocalFeaturesAndMatchesDatabase::PutFeatures(
    const std::string& image_name, const KeypointsAndDescriptors& features) {
  const std::string features_file =
      FeatureFilenameFromImage(directory_, image_name);
  CHECK(WriteKeypointsAndDescriptors(
      features_file, features.keypoints, features.descriptors))
      << "Could not write features for image " << image_name << " to file "
      << features_file;
  image_names_.insert(image_name);
  features_cache_->Insert(image_name, features);
}

// Supply an iterator to iterate over the features.
// Supply an iterator to iterate over the features.
std::vector<std::string> LocalFeaturesAndMatchesDatabase::ImageNamesOfFeatures()
    const {
  std::vector<std::string> image_names(image_names_.begin(),
                                       image_names_.end());
  return image_names;
}

size_t LocalFeaturesAndMatchesDatabase::NumImages() const {
  return image_names_.size();
}

// Get the image pair match for the images.
ImagePairMatch LocalFeaturesAndMatchesDatabase::GetImagePairMatch(
    const std::string& image_name1, const std::string& image_name2) {
  const int match_index = FindOrDieNoPrint(
      matches_index_, std::make_pair(image_name1, image_name2));
  return matches_[match_index];
}

// Set the image pair match for the images.
void LocalFeaturesAndMatchesDatabase::PutImagePairMatch(
    const std::string& image_name1,
    const std::string& image_name2,
    const ImagePairMatch& matches) {
  const auto name_pair = std::make_pair(image_name1, image_name2);

  // If the match has already been added, perform an update in place.
  if (ContainsKey(matches_index_, name_pair)) {
    matches_[matches_index_[name_pair]] = matches;
    return;
  }

  // Otherwise add the match, potentially resizing the matches container.
  std::lock_guard<std::mutex> lock(mutex_);
  matches_index_[name_pair] = matches_.size();
  matches_.push_back(matches);
}

KeypointsAndDescriptors LocalFeaturesAndMatchesDatabase::FetchImages(
    const std::string& image_name) {
  KeypointsAndDescriptors features;
  CHECK(ReadKeypointsAndDescriptors(
      FeatureFilenameFromImage(directory_, image_name),
      &features.keypoints,
      &features.descriptors));
  features.image_name = image_name;
  return features;
}

std::vector<std::pair<std::string, std::string>>
LocalFeaturesAndMatchesDatabase::ImageNamesOfMatches() const {
  std::vector<std::pair<std::string, std::string>> match_keys;
  match_keys.reserve(matches_index_.size());
  for (const auto& match : matches_index_) {
    match_keys.push_back(match.first);
  }
  return match_keys;
}

size_t LocalFeaturesAndMatchesDatabase::NumMatches() const {
  return matches_.size();
}

// Populate this database from the input matches_file, and output the view
// names and camera intrinsics.
bool LocalFeaturesAndMatchesDatabase::ReadMatchesAndGeometry(
    const std::string& matches_file,
    std::vector<std::string>* view_names,
    std::vector<CameraIntrinsicsPrior>* camera_intrinsics_prior) {
  CHECK_NOTNULL(view_names)->clear();
  CHECK_NOTNULL(camera_intrinsics_prior)->clear();

  // Return false if the file cannot be opened.
  std::ifstream matches_reader(matches_file, std::ios::in | std::ios::binary);
  if (!matches_reader.is_open()) {
    LOG(ERROR) << "Could not open the matches file: " << matches_file
               << " for reading.";
    return false;
  }

  // Make sure that Cereal is able to finish executing before returning.
  {
    cereal::PortableBinaryInputArchive input_archive(matches_reader);
    input_archive(*view_names, *camera_intrinsics_prior, matches_);
  }

  matches_index_.reserve(matches_.size());
  for (int i = 0; i < matches_.size(); i++) {
    const auto name_pair =
        std::make_pair(matches_[i].image1, matches_[i].image2);
    matches_index_[name_pair] = i;
  }

  return true;
}

// Save the matches and geometry to disk.
bool LocalFeaturesAndMatchesDatabase::SaveMatchesAndGeometry(
    const std::string& matches_file,
    const std::vector<std::string>& view_names,
    const std::vector<CameraIntrinsicsPrior>& camera_intrinsics_prior) {
  // Return false if the file cannot be opened for writing.
  std::ofstream matches_writer(matches_file, std::ios::out | std::ios::binary);
  if (!matches_writer.is_open()) {
    LOG(ERROR) << "Could not open the matches file: " << matches_file
               << " for writing.";
    return false;
  }

  // Make sure that Cereal is able to finish executing before returning.
  {
    cereal::PortableBinaryOutputArchive output_archive(matches_writer);
    output_archive(view_names, camera_intrinsics_prior, matches_);
  }

  return true;
}
}  // namespace theia
