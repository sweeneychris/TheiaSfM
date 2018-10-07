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

#ifndef THEIA_MATCHING_LOCAL_FEATURES_AND_MATCHES_DATABASE_H_
#define THEIA_MATCHING_LOCAL_FEATURES_AND_MATCHES_DATABASE_H_

#include <mutex>  // NOLINT
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "theia/matching/features_and_matches_database.h"
#include "theia/matching/image_pair_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/util/lru_cache.h"
#include "theia/util/util.h"

namespace theia {

// A simple implementation for storing features and feature matches. A local
// filesystem and cache are used to retrieve the features efficiently. The
// matches are kept in memory. This class is guaranteed to be thread safe.
class LocalFeaturesAndMatchesDatabase : public FeaturesAndMatchesDatabase {
 public:
  LocalFeaturesAndMatchesDatabase(const std::string& directory,
                                  const int max_cache_entries);
  ~LocalFeaturesAndMatchesDatabase();

  bool ContainsFeatures(const std::string& image_name) const override;

  // Get/set the features for the image. Returns true if the features exist in
  // the database and false otherwise.
  KeypointsAndDescriptors GetFeatures(const std::string& image_name) override;

  // Set the features for the image.
  void PutFeatures(const std::string& image_name,
                   const KeypointsAndDescriptors& features) override;

  // Supply an iterator to iterate over the features.
  std::vector<std::string> ImageNamesOfFeatures() const override;
  size_t NumImages() const override;

  // Get the image pair match for the images.Returns true if the features exist
  // in the database and false otherwise.
  ImagePairMatch GetImagePairMatch(const std::string& image_name1,
                                   const std::string& image_name2) override;

  // Set the image pair match for the images.
  void PutImagePairMatch(const std::string& image_name1,
                         const std::string& image_name2,
                         const ImagePairMatch& matches) override;

  std::vector<std::pair<std::string, std::string>> ImageNamesOfMatches()
      const override;
  size_t NumMatches() const override;

  // Populate this database from the input matches_file, and output the view
  // names and camera intrinsics.
  bool ReadMatchesAndGeometry(
      const std::string& matches_file,
      std::vector<std::string>* view_names,
      std::vector<CameraIntrinsicsPrior>* camera_intrinsics_prior) override;

  // Save the matches and geometry to disk.
  bool SaveMatchesAndGeometry(const std::string& matches_file,
                              const std::vector<std::string>& view_names,
                              const std::vector<CameraIntrinsicsPrior>&
                                  camera_intrinsics_prior) override;

 private:
  DISALLOW_COPY_AND_ASSIGN(LocalFeaturesAndMatchesDatabase);
  using LRUFeatureCache = LRUCache<std::string, KeypointsAndDescriptors>;

  KeypointsAndDescriptors FetchImages(const std::string& image_name);

  std::string directory_;
  std::unique_ptr<LRUFeatureCache> features_cache_;
  std::unordered_set<std::string> image_names_;
  std::unordered_map<std::pair<std::string, std::string>, int> matches_index_;
  std::vector<ImagePairMatch> matches_;
  std::mutex mutex_;
};
}  // namespace theia
#endif  // THEIA_MATCHING_LOCAL_FEATURES_AND_MATCHES_DATABASE_H_
