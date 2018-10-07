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

#ifndef THEIA_MATCHING_FEATURES_AND_MATCHES_DATABASE_H_
#define THEIA_MATCHING_FEATURES_AND_MATCHES_DATABASE_H_

#include <string>
#include <utility>
#include <vector>

#include "theia/matching/image_pair_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {

// An interface for retreiving feature and match related data. This data is
// typically memory intensive so caches or database systems may be used to
// access the data more efficiently. This class is guaranteed to be thread safe.
class FeaturesAndMatchesDatabase {
 public:
  virtual ~FeaturesAndMatchesDatabase() {}

  virtual bool ContainsFeatures(const std::string& image_name) const = 0;

  // Get/set the features for the image.
  virtual KeypointsAndDescriptors GetFeatures(
      const std::string& image_name) = 0;

  // Set the features for the image.
  virtual void PutFeatures(const std::string& image_name,
                           const KeypointsAndDescriptors& features) = 0;

  // Supply an iterator to iterate over the features.
  virtual std::vector<std::string> ImageNamesOfFeatures() const = 0;
  virtual size_t NumImages() const = 0;

  // Get the image pair match for the images.
  virtual ImagePairMatch GetImagePairMatch(const std::string& image_name1,
                                           const std::string& image_name2) = 0;

  // Set the image pair match for the images.
  virtual void PutImagePairMatch(const std::string& image_name1,
                                 const std::string& image_name2,
                                 const ImagePairMatch& matches) = 0;

  // Supply an iterator to iterate over the matches.
  virtual std::vector<std::pair<std::string, std::string>> ImageNamesOfMatches()
      const = 0;
  virtual size_t NumMatches() const = 0;

  // Populate this database from the input matches_file, and output the view
  // names and camera intrinsics.
  virtual bool ReadMatchesAndGeometry(
      const std::string& matches_file,
      std::vector<std::string>* view_names,
      std::vector<CameraIntrinsicsPrior>* camera_intrinsics_prior) = 0;

  // Save the matches and geometry to disk.
  virtual bool SaveMatchesAndGeometry(
      const std::string& matches_file,
      const std::vector<std::string>& view_names,
      const std::vector<CameraIntrinsicsPrior>& camera_intrinsics_prior) = 0;
};
}  // namespace theia
#endif  // THEIA_MATCHING_FEATURES_AND_MATCHES_DATABASE_H_
