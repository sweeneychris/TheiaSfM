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

#ifndef THEIA_SFM_RECONSTRUCTION_BUILDER_H_
#define THEIA_SFM_RECONSTRUCTION_BUILDER_H_

#include <memory>
#include <string>
#include <vector>

#include "theia/io/write_matches.h"
#include "theia/util/util.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/exif_reader.h"
#include "theia/sfm/feature_extractor.h"
#include "theia/sfm/match_and_verify_features.h"
#include "theia/sfm/reconstruction_estimator_options.h"

namespace theia {

class Reconstruction;
class TrackBuilder;
class ViewGraph;

struct ReconstructionBuilderOptions {
  // Number of threads used. Each stage of the pipeline (feature extraction,
  // matching, estimation, etc.) will use this number of threads.
  int num_threads = 1;

  // By default, the ReconstructionBuilder will attempt to reconstruct as many
  // models as possible from the input data. If set to true, only the largest
  // connected component is reconstructed.
  bool reconstruct_largest_connected_component = false;

  // Set to true to only accept calibrated views (from EXIF or elsewhere) as
  // valid inputs to the reconstruction process. When uncalibrated views are
  // added to the reconstruction builder they are ignored with a LOG warning.
  bool only_calibrated_views = false;

  // Maximum allowable track length. Tracks that are too long are exceedingly
  // likely to contain outliers.
  int max_track_length = 20;

  // Minimum number of geometrically verified inliers that a view pair must have
  // in order to be considered a good match.
  int min_num_inlier_matches = 30;

  // Descriptor type for extracting features.
  // See //theia/sfm/feature_extractor.h
  DescriptorExtractorType descriptor_type = DescriptorExtractorType::SIFT;

  // Matching strategy type.
  // See //theia/sfm/match_and_verify_features.h
  MatchingStrategy matching_strategy = MatchingStrategy::BRUTE_FORCE;

  // Options for computing matches between images.
  // See //theia/matching/feature_matcher_options.h
  FeatureMatcherOptions matching_options;

  // Settings for estimating the relative pose between two images to perform
  // geometric verification.
  // See //theia/sfm/verify_two_view_matches.h
  VerifyTwoViewMatchesOptions geometric_verification_options;

  // Options for estimating the reconstruction.
  // See //theia/sfm/reconstruction_estimator_options.h
  ReconstructionEstimatorOptions reconstruction_estimator_options;

  // If you want the matches to be saved, set this variable to the filename that
  // you want the matches to be written to. Image names, inlier matches, and
  // view metadata so that the view graph and tracks may be exactly
  // recreated.
  std::string output_matches_file;
};

// Base class for building SfM reconstructions. This class will manage the
// entire reconstruction estimation process.
class ReconstructionBuilder {
 public:
  explicit ReconstructionBuilder(const ReconstructionBuilderOptions& options);
  ~ReconstructionBuilder();

  // Add an image to the reconstruction.
  bool AddImage(const std::string& image_filepath);

  // Same as above, but with the camera priors manually specified.
  bool AddImageWithCameraIntrinsicsPrior(
      const std::string& image_filepath,
      const CameraIntrinsicsPrior& camera_intrinsics_prior);

  // Add a match to the view graph. Either this method is repeatedly called or
  // ExtractAndMatchFeatures must be called.
  bool AddTwoViewMatch(const std::string& image1,
                       const std::string& image2,
                       const ImagePairMatch& matches);

  // Extracts features and performs matching with geometric verification.
  bool ExtractAndMatchFeatures();

  // Initializes the reconstruction and view graph explicitly. This method
  // should be used as an alternative to the Add* methods.
  //
  // NOTE: The ReconstructionBuilder takses ownership of the reconstruction and
  // view graph.
  void InitializeReconstructionAndViewGraph(Reconstruction* reconstruction,
                                            ViewGraph* view_graph);

  // Estimates a Structure-from-Motion reconstruction using the specified
  // ReconstructionEstimator. Features are first extracted and matched if
  // necessary, then a reconstruction is estimated. Once a reconstruction has
  // been estimated, all views that have been successfully estimated are added
  // to the output vector and we estimate a reconstruction from the remaining
  // unestimated views. We repeat this process until no more views can be
  // successfully estimated.
  bool BuildReconstruction(std::vector<Reconstruction*>* reconstructions);

 private:
  // Matches features and adds the matches to the view graph and feature
  // correspondences to the track builder.
  template <class DescriptorType>
  bool MatchFeatures(
      const std::vector<std::vector<Keypoint> >& keypoints,
      const std::vector<std::vector<DescriptorType> >& descriptors);

  // Adds the given matches as edges in the view graph.
  void AddMatchToViewGraph(const ViewId view_id1,
                           const ViewId view_id2,
                           const ImagePairMatch& image_matches);

  // Builds tracks from the two view inlier correspondences after geometric
  // verification.
  void AddTracksForMatch(const ViewId view_id1,
                         const ViewId view_id2,
                         const ImagePairMatch& image_matches);

  // Removes all uncalibrated views from the reconstruction and view graph.
  void RemoveUncalibratedViews();

  const ReconstructionBuilderOptions options_;

  // SfM objects.
  std::unique_ptr<TrackBuilder> track_builder_;
  std::unique_ptr<Reconstruction> reconstruction_;
  std::unique_ptr<ViewGraph> view_graph_;

  // Container of image information.
  std::vector<std::string> image_filepaths_;
  std::vector<CameraIntrinsicsPrior> camera_intrinsics_priors_;

  // Exif reader.
  std::unique_ptr<ExifReader> exif_reader_;

  DISALLOW_COPY_AND_ASSIGN(ReconstructionBuilder);
};

template <class DescriptorType>
bool ReconstructionBuilder::MatchFeatures(
    const std::vector<std::vector<Keypoint> >& keypoints,
    const std::vector<std::vector<DescriptorType> >& descriptors) {
  CHECK_EQ(image_filepaths_.size(), keypoints.size());
  CHECK_EQ(descriptors.size(), keypoints.size());

  // Set up options.
  MatchAndVerifyFeaturesOptions match_and_verify_features_options;
  match_and_verify_features_options.num_threads = options_.num_threads;
  match_and_verify_features_options.min_num_inlier_matches =
      options_.min_num_inlier_matches;
  match_and_verify_features_options.matching_strategy =
      options_.matching_strategy;
  match_and_verify_features_options.feature_matcher_options =
      options_.matching_options;
  match_and_verify_features_options.geometric_verification_options =
      options_.geometric_verification_options;
  match_and_verify_features_options.geometric_verification_options
      .min_num_inlier_matches = options_.min_num_inlier_matches;

  // Match images and perform geometric verification.
  std::vector<ImagePairMatch> matches;
  CHECK(MatchAndVerifyFeatures(match_and_verify_features_options,
                               camera_intrinsics_priors_,
                               keypoints,
                               descriptors,
                               &matches)) << "Could not match features.";
  const int num_total_view_pairs =
      keypoints.size() * (keypoints.size() - 1) / 2;
  LOG(INFO) << matches.size() << " of " << num_total_view_pairs
            << " view pairs were matched and geometrically verified.";

  // Get the image names to use as the unique view name.
  std::vector<std::string> image_filenames(image_filepaths_.size());
  for (int i = 0; i < image_filenames.size(); i++) {
    const std::string& image_filepath = image_filepaths_[i];
    CHECK(GetFilenameFromFilepath(image_filepath, true, &image_filenames[i]));
  }

  // Write the matches to a file if it exists.
  if (options_.output_matches_file.length() > 0) {
    LOG(INFO) << "Writing matches to file: " << options_.output_matches_file;
    CHECK(WriteMatchesAndGeometry(options_.output_matches_file,
                                  image_filenames,
                                  camera_intrinsics_priors_,
                                  matches))
        << "Could not write the matches to " << options_.output_matches_file;
  }

  // Add the matches to the view graph and reconstruction.
  for (const auto& match : matches) {
    AddTwoViewMatch(image_filenames[match.image1_index],
                    image_filenames[match.image2_index],
                    match);
  }

  return true;
}

}  // namespace theia

#endif  // THEIA_SFM_RECONSTRUCTION_BUILDER_H_
