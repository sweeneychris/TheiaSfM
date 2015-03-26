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

#ifndef THEIA_SFM_FEATURE_EXTRACTOR_H_
#define THEIA_SFM_FEATURE_EXTRACTOR_H_

#include <Eigen/Core>
#include <algorithm>
#include <string>
#include <vector>

#include "theia/alignment/alignment.h"
#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/image.h"
#include "theia/image/keypoint_detector/sift_parameters.h"
#include "theia/util/filesystem.h"
#include "theia/util/threadpool.h"
#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {

// The various types of feature descriptors you can choose. We use the default
// keypoint extractor for each feature type. Since this is a convenience class
// anyways, this functionality is acceptable. If more flexibility (custom
// features and custom descriptors) is needed then a new class may be developed.
enum class DescriptorExtractorType {
  SIFT = 0,
  BRIEF = 1,
  BRISK = 2,
  FREAK = 3
};

struct FeatureExtractorOptions {
  int num_threads = 1;
  DescriptorExtractorType descriptor_extractor_type =
      DescriptorExtractorType::SIFT;
  // Sift parameters.
  SiftParameters sift_parameters;
  // The features returned will be no larger than this size.
  int max_num_features = 16384;
};

// Reads in the set of images provided then extracts descriptors using the
// desired descriptor type. This method can be run with multiple threads.
class FeatureExtractor {
 public:
  explicit FeatureExtractor(const FeatureExtractorOptions& options)
      : options_(options) {}
  ~FeatureExtractor() {}

  // Method to extract descriptors. Descriptors must be a float descriptor
  // Eigen::VectorXf (e.g., SIFT) or a BinaryVectorX.
  template <typename DescriptorType>
  bool Extract(const std::vector<std::string>& filenames,
               std::vector<std::vector<Keypoint> >* keypoints,
               std::vector<std::vector<DescriptorType> >* descriptors);

 private:
  // Extracts the features and metadata for a single image. This function is
  // called by the threadpool and is thus thread safe.
  template <typename DescriptorType>
  bool ExtractFeatures(const std::string& filename,
                       std::vector<Keypoint>* keypoints,
                       std::vector<DescriptorType>* descriptors);

  // Factory method to create the keypoint detector and descriptor extractor.
  std::unique_ptr<DescriptorExtractor> CreateDescriptorExtractor(
      const DescriptorExtractorType& descriptor_extractor_type);

  const FeatureExtractorOptions options_;

  DISALLOW_COPY_AND_ASSIGN(FeatureExtractor);
};

// ---------------------- Implementation ------------------------- //
template <typename DescriptorType>
bool FeatureExtractor::Extract(
    const std::vector<std::string>& filenames,
    std::vector<std::vector<Keypoint> >* keypoints,
    std::vector<std::vector<DescriptorType> >* descriptors) {
  CHECK_NOTNULL(keypoints)->resize(filenames.size());
  CHECK_NOTNULL(descriptors)->resize(filenames.size());

  // The thread pool will wait to finish all jobs when it goes out of scope.
  const int num_threads =
      std::min(options_.num_threads, static_cast<int>(filenames.size()));
  ThreadPool feature_extractor_pool(num_threads);
  for (int i = 0; i < filenames.size(); i++) {
    if (!FileExists(filenames[i])) {
      LOG(ERROR) << "Could not extract features for " << filenames[i]
                 << " because the file cannot be found.";
      continue;
    }

    feature_extractor_pool.Add(
        &FeatureExtractor::ExtractFeatures<DescriptorType>,
        this,
        filenames[i],
        &(*keypoints)[i],
        &(*descriptors)[i]);
  }
  return true;
}

template <typename DescriptorType>
bool FeatureExtractor::ExtractFeatures(
    const std::string& filename,
    std::vector<Keypoint>* keypoints,
    std::vector<DescriptorType>* descriptors) {
  std::unique_ptr<FloatImage> image(new FloatImage(filename));

  // We create these variable here instead of upon the construction of the
  // object so that they can be thread-safe. We *should* be able to use the
  // static thread_local keywords, but apparently Mac OS-X's version of clang
  // does not actually support it!
  //
  // TODO(cmsweeney): Change this so that each thread in the threadpool receives
  // exactly one object.
  std::unique_ptr<DescriptorExtractor> descriptor_extractor =
      CreateDescriptorExtractor(options_.descriptor_extractor_type);

  // Exit if the descriptor extraction fails.
  if (!descriptor_extractor->DetectAndExtractDescriptors(*image,
                                                         keypoints,
                                                         descriptors)) {
    LOG(ERROR) << "Could not extract descriptors in image " << filename;
    return false;
  }

  if (keypoints->size() > options_.max_num_features) {
    keypoints->resize(options_.max_num_features);
    descriptors->resize(options_.max_num_features);
  }

  VLOG(1) << "Successfully extracted " << descriptors->size()
          << " features from image " << filename;
  return true;
}


}  // namespace theia

#endif  // THEIA_SFM_FEATURE_EXTRACTOR_H_
