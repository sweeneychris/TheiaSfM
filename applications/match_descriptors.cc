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
    descriptor, "SIFT",
    "Type of feature descriptor to use. Must be one of the following: "
    "SIFT");
DEFINE_string(matcher, "brute_force",
              "Matching type to use. Must be brute_force, or "
              "cascade_hashing");
DEFINE_double(lowes_ratio, 0.8, "Lowes ratio used to filter feature matching.");
DEFINE_int32(num_threads, 1,
             "Number of threads to use for feature extraction and matching.");
DEFINE_string(img_output_dir, ".", "Name of output image file.");

theia::DescriptorExtractorType GetDescriptorExtractorType(
    const std::string& descriptor) {
  if (descriptor == "SIFT") {
    return theia::DescriptorExtractorType::SIFT;
  } else {
    LOG(ERROR) << "Invalid DescriptorExtractor specified. Using SIFT instead.";
    return theia::DescriptorExtractorType::SIFT;
  }
}

template <class DistanceMetric>
theia::FeatureMatcher<DistanceMetric>* CreateMatcher(
    const std::string& matcher) {
  if (matcher == "cascade_hashing") {
    return new theia::CascadeHashingFeatureMatcher;
  } else if (matcher == "brute_force") {
    return new theia::BruteForceFeatureMatcher<DistanceMetric>;
  }

  LOG(ERROR) << "Invalid matcher specified";
  return nullptr;
}

template <class DistanceMetric, class DescriptorType>
void ExtractAndMatchFeatures(
    std::vector<theia::FloatImage*>* images,
    std::vector<std::vector<theia::Keypoint> >* keypoints,
    std::vector<theia::ImagePairMatch>* matches) {
  // Get image filenames and read in images.
  std::vector<std::string> img_filepaths;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_input_imgs, &img_filepaths));
  std::vector<std::vector<DescriptorType> > descriptors;
  for (int i = 0; i < img_filepaths.size(); i++) {
    images->emplace_back(new theia::FloatImage(img_filepaths[i]));
  }

  // Extract features.
  theia::FeatureExtractorOptions feature_extractor_options;
  feature_extractor_options.descriptor_extractor_type =
      GetDescriptorExtractorType(FLAGS_descriptor);
  feature_extractor_options.num_threads = FLAGS_num_threads;
  theia::FeatureExtractor feature_extractor(feature_extractor_options);
  auto start = std::chrono::system_clock::now();
  CHECK(feature_extractor.Extract(img_filepaths, keypoints, &descriptors))
      << "Feature extraction failed!";
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start);
  const double time_to_extract_features = duration.count();

  // Create feature matcher and set up options.
  std::unique_ptr<theia::FeatureMatcher<DistanceMetric> > matcher(
      CHECK_NOTNULL(CreateMatcher<DistanceMetric>(FLAGS_matcher)));
  theia::FeatureMatcherOptions options;
  options.num_threads = FLAGS_num_threads;
  options.lowes_ratio = FLAGS_lowes_ratio;

  // Match features.
  start = std::chrono::system_clock::now();
  for (int i = 0; i < keypoints->size(); i++) {
    matcher->AddImage(&keypoints->at(i), &descriptors[i]);
  }
  theia::VerifyTwoViewMatchesOptions verification_options;
  matcher->MatchImagesWithGeometricVerification(options,
                                                verification_options,
                                                matches);
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start);
  const double time_to_match_features = duration.count();

  LOG(INFO) << "It took " << (time_to_extract_features / 1000.0)
            << " seconds to extract features and "
            << (time_to_match_features / 1000.0) << " seconds to match "
            << matches->size() << " image pairs";
}

int main(int argc, char *argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  std::vector<theia::FloatImage*> images;
  std::vector<std::vector<theia::Keypoint>> keypoints;
  std::vector<theia::ImagePairMatch> image_pair_matches;
  ExtractAndMatchFeatures<theia::L2, Eigen::VectorXf>(&images,
                                                      &keypoints,
                                                      &image_pair_matches);
  for (int i = 0; i < image_pair_matches.size(); i++) {
    theia::ImageCanvas image_canvas;
    const int img1_index = image_pair_matches[i].image1_index;
    const int img2_index = image_pair_matches[i].image2_index;
    image_canvas.AddImage(*images[img1_index]);
    image_canvas.AddImage(*images[img2_index]);
    const std::string match_output = theia::StringPrintf(
        "%s/matches_%i_%i.png",
        FLAGS_img_output_dir.c_str(),
        img1_index,
        img2_index);
    image_canvas.DrawMatchedFeatures(0,
                                     1,
                                     image_pair_matches[i].correspondences,
                                     0.1);
    image_canvas.Write(match_output);
  }

  theia::STLDeleteElements(&images);
}
