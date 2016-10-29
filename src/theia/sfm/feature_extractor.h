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
#include <string>

#include "theia/alignment/alignment.h"
#include "theia/image/descriptor/create_descriptor_extractor.h"
#include "theia/util/util.h"

namespace theia {
class Keypoint;

// Reads in the set of images provided then extracts descriptors using the
// desired descriptor type. This method can be run with multiple threads.
class FeatureExtractor {
 public:
  struct Options {
    int num_threads = 1;
    DescriptorExtractorType descriptor_extractor_type =
        DescriptorExtractorType::SIFT;

    // The density of features to extract. DENSE means more features are
    // extracted per image and SPARSE means fewer features per image are
    // extracted.
    FeatureDensity feature_density = FeatureDensity::NORMAL;

    // The features returned will be no larger than this size.
    int max_num_features = 16384;

    // If we wish to write the features to disk, they will be output in this
    // directory with the same name as the input image and a ".features"
    // appended.
    std::string output_directory = "";
  };

  explicit FeatureExtractor(const Options& options)
      : options_(options), write_features_to_disk_(false) {}
  ~FeatureExtractor() {}

  // Method to extract descriptors.
  bool Extract(const std::vector<std::string>& filenames,
               std::vector<std::vector<Keypoint> >* keypoints,
               std::vector<std::vector<Eigen::VectorXf> >* descriptors);

  // Extracts descriptors and writes them to disk. The features from each image
  // are written to individual files in the directory specified in the options.
  bool ExtractToDisk(const std::vector<std::string>& filenames);

 private:
  // Extracts the features and metadata for a single image. This function is
  // called by the threadpool and is thus thread safe.
  bool ExtractFeatures(const std::string& filename,
                       std::vector<Keypoint>* keypoints,
                       std::vector<Eigen::VectorXf>* descriptors);

  const Options options_;
  bool write_features_to_disk_;

  DISALLOW_COPY_AND_ASSIGN(FeatureExtractor);
};

}  // namespace theia

#endif  // THEIA_SFM_FEATURE_EXTRACTOR_H_
