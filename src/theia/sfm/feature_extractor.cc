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

#include "theia/sfm/feature_extractor.h"

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/descriptor/sift_descriptor.h"
#include "theia/image/image.h"
#include "theia/util/filesystem.h"
#include "theia/util/threadpool.h"

namespace theia {

bool FeatureExtractor::Extract(
    const std::vector<std::string>& filenames,
    std::vector<std::vector<Keypoint> >* keypoints,
    std::vector<std::vector<Eigen::VectorXf> >* descriptors) {
  CHECK_GT(filenames.size(), 0) << "FeatureExtractor::Extract requires at "
                                   "least one image in order to extract "
                                   "features.";
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
        &FeatureExtractor::ExtractFeatures,
        this,
        filenames[i],
        &(*keypoints)[i],
        &(*descriptors)[i]);
  }
  return true;
}

bool FeatureExtractor::ExtractFeatures(
    const std::string& filename,
    std::vector<Keypoint>* keypoints,
    std::vector<Eigen::VectorXf>* descriptors) {
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

std::unique_ptr<DescriptorExtractor>
FeatureExtractor::CreateDescriptorExtractor(
    const DescriptorExtractorType& descriptor_extractor_type) {
  std::unique_ptr<DescriptorExtractor> descriptor_extractor;
  switch (descriptor_extractor_type) {
    case DescriptorExtractorType::SIFT:
      descriptor_extractor.reset(
          new SiftDescriptorExtractor(options_.sift_parameters));
      break;
    default:
      LOG(ERROR) << "Invalid Descriptor Extractor specified.";
  }
  CHECK(descriptor_extractor->Initialize())
      << "Could not initialize the Descriptor Extractor";
  return descriptor_extractor;
}

}  // namespace theia
