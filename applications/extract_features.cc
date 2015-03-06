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
#include <time.h>
#include <theia/theia.h>
#include <chrono>
#include <string>
#include <vector>

DEFINE_string(
    input_imgs, "",
    "Filepath of the images you want to extract features and compute matches "
    "for. The filepath should be a wildcard to match multiple images.");
DEFINE_string(
    img_output_dir, ".",
    "Name of output image directory. No trailing slash should be given.");
DEFINE_int32(num_threads, 1,
             "Number of threads to use for feature extraction and matching.");
DEFINE_string(
    descriptor, "SIFT",
    "Type of feature descriptor to use. Must be one of the following: "
    "SIFT, BRIEF, BRISK, FREAK");

theia::DescriptorExtractorType GetDescriptorExtractorType(
    const std::string& descriptor) {
  if (descriptor == "SIFT") {
    return theia::DescriptorExtractorType::SIFT;
  } else if (descriptor == "BRIEF") {
    return theia::DescriptorExtractorType::BRIEF;
  } else if (descriptor == "BRISK") {
    return theia::DescriptorExtractorType::BRISK;
  } else if (descriptor == "FREAK") {
    return theia::DescriptorExtractorType::FREAK;
  } else {
    LOG(ERROR) << "Invalid DescriptorExtractor specified. Using SIFT instead.";
    return theia::DescriptorExtractorType::SIFT;
  }
}

int main(int argc, char *argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  // Get image filenames.
  std::vector<std::string> img_filepaths;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_input_imgs, &img_filepaths));

  // Set up the feature extractor.
  theia::FeatureExtractorOptions feature_extractor_options;
  feature_extractor_options.descriptor_extractor_type =
      GetDescriptorExtractorType(FLAGS_descriptor);
  feature_extractor_options.num_threads = FLAGS_num_threads;
  theia::FeatureExtractor feature_extractor(feature_extractor_options);

  std::vector<std::vector<theia::Keypoint>* > keypoints;
  std::vector<std::vector<Eigen::VectorXf>* > descriptors;
  std::vector<std::vector<theia::BinaryVectorX>* > binary_descriptors;

  // Extract features from all images.
  double time_to_extract_features;
  if (FLAGS_descriptor == "SIFT") {
    auto start = std::chrono::system_clock::now();
    CHECK(feature_extractor.Extract(img_filepaths,
                                    &keypoints,
                                    &descriptors))
        << "Feature extraction failed!";
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start);
    time_to_extract_features = duration.count();
  } else {
    auto start = std::chrono::system_clock::now();
    CHECK(feature_extractor.Extract(img_filepaths,
                                    &keypoints,
                                    &binary_descriptors))
        << "Feature extraction failed!";
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start);
    time_to_extract_features = duration.count();
  }

  for (int i = 0; i < keypoints.size(); i++) {
    theia::FloatImage image(img_filepaths[i]);
    theia::ImageCanvas image_canvas;
    image_canvas.AddImage(image);
    const std::string feature_output = theia::StringPrintf(
        "%s/detected_image_%i.png", FLAGS_img_output_dir.c_str(), i);
    const theia::RGBPixel color = {255.0, 0.0, 0.0};
    image_canvas.DrawFeatures(*keypoints[i], color);
    image_canvas.Write(feature_output);
  }

  LOG(INFO) << "It took " << (time_to_extract_features / 1000.0)
            << " seconds to extract descriptors from " << img_filepaths.size()
            << " images.";

  theia::STLDeleteElements(&keypoints);
  theia::STLDeleteElements(&descriptors);
  theia::STLDeleteElements(&binary_descriptors);
}
