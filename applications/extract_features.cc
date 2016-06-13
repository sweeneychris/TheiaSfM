// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>
#include <string>
#include <vector>

#include "applications/command_line_helpers.h"

DEFINE_string(
    input_images, "",
    "Filepath of the images you want to extract features and compute matches "
    "for. The filepath should be a wildcard to extract features from multiple "
    "images.");
DEFINE_string(features_output_directory, ".",
              "Name of output directory to write the features files.");
DEFINE_int32(num_threads, 1,
             "Number of threads to use for feature extraction and matching.");
DEFINE_string(
    descriptor, "SIFT",
    "Type of feature descriptor to use. Must be one of the following: "
    "SIFT");
DEFINE_string(feature_density, "NORMAL",
              "Set to SPARSE, NORMAL, or DENSE to extract fewer or more "
              "features from each image.");

int main(int argc, char *argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  // Get image filenames.
  std::vector<std::string> img_filepaths;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_input_images, &img_filepaths));
  CHECK_GT(img_filepaths.size(), 0)
      << "No images found in: " << FLAGS_input_images;

  // Set up the feature extractor.
  theia::FeatureExtractor::Options options;
  options.descriptor_extractor_type =
      StringToDescriptorExtractorType(FLAGS_descriptor);
  options.feature_density = StringToFeatureDensity(FLAGS_feature_density);
  options.num_threads = FLAGS_num_threads;
  options.output_directory = FLAGS_features_output_directory;

  theia::FeatureExtractor feature_extractor(options);

  // Extract features from all images.
  theia::Timer timer;
  CHECK(feature_extractor.ExtractToDisk(img_filepaths))
      << "Feature extraction failed!";
  const double time_to_extract_features = timer.ElapsedTimeInSeconds();

  LOG(INFO) << "It took " << time_to_extract_features
            << " seconds to extract descriptors from " << img_filepaths.size()
            << " images.";
}
