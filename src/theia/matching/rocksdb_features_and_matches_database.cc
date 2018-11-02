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

#include "theia/matching/rocksdb_features_and_matches_database.h"

#include <cstdlib>
#include <glog/logging.h>
#include <istream>
#include <mutex>
#include <streambuf>
#include <string>

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <rocksdb/db.h>
#include <rocksdb/filter_policy.h>
#include <rocksdb/table.h>

#include "theia/matching/image_pair_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/util/filesystem.h"
#include "theia/util/map_util.h"
#include "theia/util/string.h"

namespace theia {
namespace {
using StringPair = std::pair<std::string, std::string>;

static const std::string kFeaturesColumnFamilyName =
    "keypoints_and_descriptors";
static const std::string kMatchesColumnFamilyName = "image_pair_matches";
static const std::string kIntrinsicsColumnFamilyName =
    "camera_intrinsics_prior";
static const std::string kNamePairSeparator = "/";

// For serialization using the Cereal library we must provide a stream for the
// data. This struct allows for the results from RocksDB to be directly consumed
// by Cereal without having to copy the data.
struct ZeroCopyBuffer : std::streambuf {
  ZeroCopyBuffer(char const* base, size_t size) {
    char* p(const_cast<char*>(base));
    this->setg(p, p, p + size);
  }
};

// Creates a column family with the specified name and returns the handle.s
rocksdb::ColumnFamilyHandle* CreateColumnFamily(const rocksdb::Options& options,
                                                const std::string& column_name,
                                                rocksdb::DB* database) {
  rocksdb::ColumnFamilyHandle* temp_col_family_handle = nullptr;
  database->CreateColumnFamily(options, column_name, &temp_col_family_handle);
  return temp_col_family_handle;
}

std::string ComposeImageNamePair(const std::string& image1,
                                 const std::string& image2) {
  return image1 + kNamePairSeparator + image2;
}

StringPair DecomposeImageNamePair(const std::string& image_pair) {
  const auto delimiter_index = image_pair.find(kNamePairSeparator);
  return std::make_pair(image_pair.substr(0, delimiter_index),
                        image_pair.substr(delimiter_index + 1));
}
}  // namespace

RocksDbFeaturesAndMatchesDatabase::RocksDbFeaturesAndMatchesDatabase(
    const std::string& directory)
    : directory_(directory) {
  AppendTrailingSlashIfNeeded(&directory_);
  InitializeRocksDB();
}

void RocksDbFeaturesAndMatchesDatabase::InitializeRocksDB() {
  options_.reset(new rocksdb::Options);
  // Number of threads for writing to disk.
  options_->max_background_jobs = 4;
  options_->db_write_buffer_size = 1 << 30;
  // 1 MB.
  options_->bytes_per_sync = 1 << 20;
  options_->compaction_pri = rocksdb::kMinOverlappingRatio;
  options_->create_if_missing = true;
  options_->level_compaction_dynamic_level_bytes = true;
  options_->statistics = rocksdb::CreateDBStatistics();

  rocksdb::BlockBasedTableOptions table_options;
  // 512 MB.
  table_options.block_cache = rocksdb::NewLRUCache(512 << 20);
  table_options.block_size = 16 * 1024;
  table_options.cache_index_and_filter_blocks = true;
  table_options.pin_l0_filter_and_index_blocks_in_cache = true;
  table_options.filter_policy.reset(rocksdb::NewBloomFilterPolicy(10, false));
  options_->table_factory.reset(
      rocksdb::NewBlockBasedTableFactory(table_options));

  // Get column family descriptors to open the database.
  std::vector<rocksdb::ColumnFamilyDescriptor> column_descriptors;

  // Get the column families.
  std::vector<std::string> existing_column_families;
  rocksdb::Status status = rocksdb::DB::ListColumnFamilies(
      *options_, directory_, &existing_column_families);
  // If the database exists, retain the column family descriptors.
  if (status.ok()) {
    LOG(INFO) << "Reading existing DB.";
    // Create column descriptors from the existing names.
    for (const std::string& column_family : existing_column_families) {
      column_descriptors.emplace_back(column_family, *options_);
    }
  } else {
    // RocksDB requires you to have the default column family when creating a
    // database.
    column_descriptors.emplace_back(rocksdb::kDefaultColumnFamilyName,
                                    *options_);
  }

  // Open the DB, creating it if necessary.
  rocksdb::DB* temp_db;
  std::vector<rocksdb::ColumnFamilyHandle*> temp_col_family_handles;
  status = rocksdb::DB::Open(*options_,
                             directory_,
                             column_descriptors,
                             &temp_col_family_handles,
                             &temp_db);
  CHECK(status.ok()) << "RocksDB could not open the database at: " << directory_
                     << "\n"
                     << status.ToString();

  // Take ownership of the database object.
  database_.reset(temp_db);

  // If the DB is new (i.e. no existing column families are present), then we
  // need to create the column families.
  if (existing_column_families.empty()) {
    features_handle_.reset(CreateColumnFamily(
        *options_, kFeaturesColumnFamilyName, database_.get()));
    matches_handle_.reset(CreateColumnFamily(
        *options_, kMatchesColumnFamilyName, database_.get()));
    intrinsics_prior_handle_.reset(CreateColumnFamily(
        *options_, kIntrinsicsColumnFamilyName, database_.get()));
  } else {
    // Otherwise, set up the mapping for the existing column families in the
    // database.
    for (int i = 0; i < temp_col_family_handles.size(); i++) {
      if (existing_column_families[i] == kFeaturesColumnFamilyName) {
        features_handle_.reset(temp_col_family_handles[i]);
      } else if (existing_column_families[i] == kMatchesColumnFamilyName) {
        matches_handle_.reset(temp_col_family_handles[i]);
      } else if (existing_column_families[i] == kIntrinsicsColumnFamilyName) {
        intrinsics_prior_handle_.reset(temp_col_family_handles[i]);
      }
    }
  }
}

RocksDbFeaturesAndMatchesDatabase::~RocksDbFeaturesAndMatchesDatabase() {}

bool RocksDbFeaturesAndMatchesDatabase::ContainsCameraIntrinsicsPrior(
    const std::string& image_name) {
  rocksdb::ReadOptions options;
  const rocksdb::Slice key(image_name);
  rocksdb::PinnableSlice value;
  const rocksdb::Status status =
      database_->Get(options, intrinsics_prior_handle_.get(), key, &value);
  return !status.IsNotFound();
}

// Get/set the features for the image.
CameraIntrinsicsPrior
RocksDbFeaturesAndMatchesDatabase::GetCameraIntrinsicsPrior(
    const std::string& image_name) {
  rocksdb::ReadOptions options;
  const rocksdb::Slice key(image_name);
  rocksdb::PinnableSlice value;
  const rocksdb::Status status =
      database_->Get(options, intrinsics_prior_handle_.get(), key, &value);
  CHECK(!status.IsNotFound())
      << "Could not find intrinsics for " << image_name << " in the database.";

  // Create a stream wrapped around the rocksdb value.
  ZeroCopyBuffer buffer(value.data(), value.size());
  std::istream ins(&buffer);

  // Load the keypoints and descriptors.
  CameraIntrinsicsPrior intrinsics_prior;
  {
    cereal::PortableBinaryInputArchive input_archive(ins);
    input_archive(intrinsics_prior);
  }

  std::string out;
  database_->GetProperty("rocksdb.estimate-table-readers-mem", &out);
  return intrinsics_prior;
}

// Set the features for the image.
void RocksDbFeaturesAndMatchesDatabase::PutCameraIntrinsicsPrior(
    const std::string& image_name, const CameraIntrinsicsPrior& intrinsics) {
  std::stringstream ss;
  {
    cereal::PortableBinaryOutputArchive output_archive(ss);
    output_archive(intrinsics);
  }
  rocksdb::WriteOptions options;
  const rocksdb::Slice key(image_name);
  const rocksdb::Status status =
      database_->Put(options, intrinsics_prior_handle_.get(), key, ss.str());
  CHECK(status.ok()) << "Could not insert intrinsics for " << image_name
                     << " into the database.";
}

// Supply an iterator to iterate over the priors.
std::vector<std::string>
RocksDbFeaturesAndMatchesDatabase::ImageNamesOfCameraIntrinsicsPriors() {
  // Iterate over the features column family and grab the keys.
  std::vector<std::string> image_names;
  auto it = database_->NewIterator(rocksdb::ReadOptions(),
                                   intrinsics_prior_handle_.get());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    image_names.push_back(it->key().ToString());
  }

  return image_names;
}

size_t RocksDbFeaturesAndMatchesDatabase::NumCameraIntrinsicsPrior() {
  std::uint64_t num_images;
  database_->GetIntProperty(
      intrinsics_prior_handle_.get(), "rocksdb.estimate-num-keys", &num_images);
  return static_cast<size_t>(num_images);
}

bool RocksDbFeaturesAndMatchesDatabase::ContainsFeatures(
    const std::string& image_name) {
  rocksdb::ReadOptions options;
  const rocksdb::Slice key(image_name);
  rocksdb::PinnableSlice value;
  const rocksdb::Status status =
      database_->Get(options, features_handle_.get(), key, &value);
  return !status.IsNotFound();
}

// Get/set the features for the image.
KeypointsAndDescriptors RocksDbFeaturesAndMatchesDatabase::GetFeatures(
    const std::string& image_name) {
  rocksdb::ReadOptions options;
  const rocksdb::Slice key(image_name);
  rocksdb::PinnableSlice value;
  const rocksdb::Status status =
      database_->Get(options, features_handle_.get(), key, &value);
  CHECK(!status.IsNotFound())
      << "Could not find features for " << image_name << " in the database.";

  // Create a stream wrapped around the rocksdb value.
  ZeroCopyBuffer buffer(value.data(), value.size());
  std::istream ins(&buffer);

  // Load the keypoints and descriptors.
  KeypointsAndDescriptors features;
  {
    cereal::PortableBinaryInputArchive input_archive(ins);
    input_archive(
        features.image_name, features.keypoints, features.descriptors);
  }
  return features;
}

// Set the features for the image.
void RocksDbFeaturesAndMatchesDatabase::PutFeatures(
    const std::string& image_name, const KeypointsAndDescriptors& features) {
  std::stringstream ss;
  {
    cereal::PortableBinaryOutputArchive output_archive(ss);
    output_archive(
        features.image_name, features.keypoints, features.descriptors);
  }

  rocksdb::WriteOptions options;
  const rocksdb::Slice key(image_name);
  const rocksdb::Status status =
      database_->Put(options, features_handle_.get(), key, ss.str());
  CHECK(status.ok()) << "Could not insert features for " << image_name
                     << " into the database.";
}

std::vector<std::string>
RocksDbFeaturesAndMatchesDatabase::ImageNamesOfFeatures() {
  // Iterate over the features column family and grab the keys.
  std::vector<std::string> image_names;
  auto it =
      database_->NewIterator(rocksdb::ReadOptions(), features_handle_.get());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    image_names.push_back(it->key().ToString());
  }

  return image_names;
}

size_t RocksDbFeaturesAndMatchesDatabase::NumImages() {
  std::uint64_t num_images;
  database_->GetIntProperty(
      features_handle_.get(), "rocksdb.estimate-num-keys", &num_images);
  return static_cast<size_t>(num_images);
}

// Get the image pair match for the images.
ImagePairMatch RocksDbFeaturesAndMatchesDatabase::GetImagePairMatch(
    const std::string& image_name1, const std::string& image_name2) {
  const std::string image_name_pair =
      ComposeImageNamePair(image_name1, image_name2);

  rocksdb::ReadOptions options;
  const rocksdb::Slice key(image_name_pair);
  rocksdb::PinnableSlice value;
  const rocksdb::Status status =
      database_->Get(options, matches_handle_.get(), key, &value);
  CHECK(!status.IsNotFound()) << "Could not find the image pair match for ("
                              << image_name1 << ", " << image_name2 << ")";

  // Create a stream wrapped around the rocksdb value.
  ZeroCopyBuffer buffer(value.data(), value.size());
  std::istream ins(&buffer);

  // Load the keypoints and descriptors.
  ImagePairMatch matches;
  {
    cereal::PortableBinaryInputArchive input_archive(ins);
    input_archive(matches);
  }
  return matches;
}

// Set the image pair match for the images.
void RocksDbFeaturesAndMatchesDatabase::PutImagePairMatch(
    const std::string& image_name1,
    const std::string& image_name2,
    const ImagePairMatch& matches) {
  const std::string image_name_pair =
      ComposeImageNamePair(image_name1, image_name2);

  std::stringstream ss;
  {
    cereal::PortableBinaryOutputArchive output_archive(ss);
    output_archive(matches);
  }

  rocksdb::WriteOptions options;
  const rocksdb::Slice key(image_name_pair);
  const rocksdb::Status status =
      database_->Put(options, matches_handle_.get(), key, ss.str());
  CHECK(status.ok());
}

std::vector<StringPair>
RocksDbFeaturesAndMatchesDatabase::ImageNamesOfMatches() {
  // Iterate over the features column family and grab the keys.
  std::vector<StringPair> image_match_names;
  auto it =
      database_->NewIterator(rocksdb::ReadOptions(), matches_handle_.get());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    image_match_names.push_back(DecomposeImageNamePair(it->key().ToString()));
  }

  return image_match_names;
}

size_t RocksDbFeaturesAndMatchesDatabase::NumMatches() {
  std::uint64_t num_matches;
  database_->GetIntProperty(
      matches_handle_.get(), "rocksdb.estimate-num-keys", &num_matches);
  return static_cast<size_t>(num_matches);
}

void RocksDbFeaturesAndMatchesDatabase::RemoveAllMatches() {
  // Drop the column family handle -- this deletes all key/values in the column
  // family.
  database_->DropColumnFamily(matches_handle_.get());

  // Add the column family back again.
  matches_handle_.reset(
      CreateColumnFamily(*options_, kMatchesColumnFamilyName, database_.get()));
}

}  // namespace theia
