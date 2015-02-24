// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#include <theia/theia.h>
#include <glog/logging.h>
#include <gflags/gflags.h>

#include <fstream>  // NOLINT
#include <string>
#include <utility>
#include <unordered_map>
#include <unordered_set>

DEFINE_string(list_file, "", "list.txt filepath");
DEFINE_string(eg_file, "", "EGs.txt filepath");
DEFINE_string(cc_file, "", "cc.txt filepath");
DEFINE_string(coords_file, "", "coords.txt filepath");
DEFINE_string(tracks_file, "", "tracks.txt filepath");
DEFINE_string(output_match_file, "", "Location of match file to write.");

// Geometric verification options.
DEFINE_int32(min_num_verified_matches, 30,
            "Minimum number of matches to be considered an inlier.");

using theia::CameraIntrinsicsPrior;
using theia::ImagePairMatch;

typedef std::pair<int, int> IntPair;

// Reads the list file and ignores the focal length since it can be recovered
// from the coords file.
bool ReadListsFile(const std::string& list_filename,
                   std::vector<std::string>* images) {
  std::ifstream ifs(list_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the list file from " << list_filename;
    return false;
  }

  while (!ifs.eof()) {
    // Read in the filename.
    std::string filename, truncated_filename;
    CameraIntrinsicsPrior intrinsics;

    ifs >> filename;
    if (filename.length() == 0) {
      break;
    }
    CHECK(theia::GetFilenameFromFilepath(filename, false, &truncated_filename));
    images->emplace_back(truncated_filename);

    // Ignore the rest of the line.
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return true;
}

// Reads the epipolar geometry files.
bool ReadEGs(const std::string& eg_filename,
             std::vector<ImagePairMatch>* matches) {
  std::ifstream ifs(eg_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the EG file from " << eg_filename;
    return false;
  }

  while (!ifs.eof()) {
    ImagePairMatch image_pair_match;

    ifs >> image_pair_match.image1_index;
    ifs >> image_pair_match.image2_index;
    CHECK_LT(image_pair_match.image1_index, image_pair_match.image2_index);

    // The rotation defines the camera 2 to camera 1 transformation in row-major
    // order). We want a camera 1 to camera 2 transformation so we read in the
    // transpose (i.e., column-major order).
    Eigen::Matrix3d rotation;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        ifs >> rotation(i, j);
      }
    }
    const Eigen::AngleAxisd rotation_aa(rotation);
    image_pair_match.twoview_info.rotation_2 =
        rotation_aa.angle() * rotation_aa.axis();

    // Assign the position.
    ifs >> image_pair_match.twoview_info.position_2[0];
    ifs >> image_pair_match.twoview_info.position_2[1];
    ifs >> image_pair_match.twoview_info.position_2[2];

    // Add the match to the output.
    matches->emplace_back(image_pair_match);
  }
  return true;
}

// Reads the connected components file.
bool ReadCC(const std::string& cc_filename,
            std::unordered_set<int>* valid_images) {
  std::ifstream ifs(cc_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the cc file from " << cc_filename;
    return false;
  }

  while (!ifs.eof()) {
    int img_index;
    ifs >> img_index;
    valid_images->insert(img_index);
  }

  return true;
}

// Coords is a map of a pair <img id, feature id> to feature coordinate.
// Tracks is a map of track_id to a vector of <img id, feature id>.
bool ReadTracks(
    const std::string& tracks_filename,
    std::unordered_map<int, std::unordered_map<int, Eigen::Vector2d> >* coords,
    std::unordered_map<int, std::vector<IntPair> >* tracks) {
  std::ifstream ifs(tracks_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the coords file from " << tracks_filename;
    return false;
  }

  // Read number of tracks.
  int num_tracks;
  ifs >> num_tracks;

  tracks->reserve(num_tracks);
  coords->reserve(num_tracks * 15);
  for (int i = 0; i < num_tracks; i++) {
    std::vector<IntPair> track;
    int num_features;
    ifs >> num_features;
    track.resize(num_features);
    for (int j = 0; j < num_features; j++) {
      ifs >> track[j].first;
      ifs >> track[j].second;

      // Place an empty feature in the coords list. This will be used later to
      // verify that the feature exists in a track. Only features that exist in
      // a track will be read from the coords file.
      (*coords)[track[j].first][track[j].second] = Eigen::Vector2d::Zero();
    }
    tracks->emplace(i, track);
  }

  return true;
}

void ReadCoordsHeaderLine(const std::string& line,
                          int* image_index,
                          int* num_keys,
                          CameraIntrinsicsPrior* prior) {
  float principal_point_x, principal_point_y, focal_length;
  char name[256];
  sscanf(line.c_str(),
         "#index = %d, name = %s keys = %d, px = %f, py = %f, focal = %f",
         image_index,
         name,
         num_keys,
         &principal_point_x,
         &principal_point_y,
         &focal_length);

  LOG(INFO) << "Image: " << name << " calibration: " << focal_length
            << " px, py = " << principal_point_x << ", " << principal_point_y;

  prior->image_width = principal_point_x * 2.0;
  prior->image_width = principal_point_y * 2.0;
  prior->focal_length.value = focal_length;
  prior->focal_length.is_set = true;
  prior->principal_point[0].value = principal_point_x;
  prior->principal_point[0].is_set = true;
  prior->principal_point[1].value = principal_point_y;
  prior->principal_point[1].is_set = true;
}

// Reads the coords file. Only the coords with a valid track in the connected
// component are kept.
bool ReadCoords(const std::string& coords_filename,
                const std::unordered_set<int>& valid_images,
                std::vector<CameraIntrinsicsPrior>* intrinsics,
                std::unordered_map<
                    int, std::unordered_map<int, Eigen::Vector2d> >* coords) {
  std::ifstream ifs(coords_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the coords file from " << coords_filename;
    return false;
  }

  while (!ifs.eof()) {
    std::string line;
    std::getline(ifs, line);
    if (ifs.eof()) {
      break;
    }

    CameraIntrinsicsPrior prior;

    int image_index, num_keys;
    ReadCoordsHeaderLine(line, &image_index, &num_keys, &prior);
    intrinsics->at(image_index) = prior;

    // Only read images with a valid connected component.
    if (!theia::ContainsKey(valid_images, image_index)) {
      for (int i = 0; i < num_keys; i++) {
          std::getline(ifs, line);
      }
      continue;
    }

    auto& feature_map = (*coords)[image_index];

    Eigen::Vector2d keypoint;
    for (int i = 0; i < num_keys; i++) {
      if (!theia::ContainsKey(feature_map, i)) {
        std::getline(ifs, line);
        continue;
      }
      std::getline(ifs, line);
      sscanf(line.c_str(), "%*d %lf %lf", &keypoint[0], &keypoint[1]);
      feature_map[i] = keypoint;
    }
  }

  return true;
}

void ConvertTracksToFeatureCorrespondences(
    const std::unordered_map<int, std::unordered_map<int, Eigen::Vector2d> >&
        coords,
    const std::unordered_map<int, std::vector<IntPair> >& tracks,
    std::vector<ImagePairMatch>* matches) {
  // Collect all the feature matches by the image id pairs.
  std::unordered_map<IntPair, std::vector<IntPair> > view_pair_to_features;
  for (const auto& track : tracks) {
    std::vector<IntPair> sorted_track = track.second;
    std::sort(sorted_track.begin(), sorted_track.end());
    // Add the minimal number of tracks.
    for (int i = 0; i < sorted_track.size() - 1; i++) {
      const IntPair image_ids(sorted_track[i].first, sorted_track[i + 1].first);
      const IntPair feature_ids(sorted_track[i].second,
                                sorted_track[i + 1].second);
      view_pair_to_features[image_ids].emplace_back(feature_ids);
    }
  }

  for (auto& match : *matches) {
    const IntPair view_id_pair(match.image1_index, match.image2_index);
    const auto* feature_matches =
        theia::FindOrNull(view_pair_to_features, view_id_pair);
    if (feature_matches == nullptr) {
      continue;
    }

    const auto& features1 = theia::FindOrDie(coords, match.image1_index);
    const auto& features2 = theia::FindOrDie(coords, match.image2_index);
    for (const auto& feature_match : *feature_matches) {
      theia::FeatureCorrespondence correspondence;
      correspondence.feature1 =
          theia::FindOrDie(features1, feature_match.first);
      correspondence.feature2 =
          theia::FindOrDie(features2, feature_match.second);
      match.correspondences.emplace_back(correspondence);
    }
  }
}

void RemoveUncalibratedImageAndReindex(
    const std::unordered_set<int>& valid_images,
    std::vector<std::string>* image_files,
    std::vector<CameraIntrinsicsPrior>* intrinsics,
    std::vector<ImagePairMatch>* matches) {
  std::vector<std::string> calibrated_image_files(valid_images.size());
  std::vector<CameraIntrinsicsPrior> calibrated_intrinsics(valid_images.size());
  std::unordered_map<int, int> old_to_new_index;
  int new_index = 0;
  for (const int old_index : valid_images) {
    old_to_new_index[old_index] = new_index;
    calibrated_image_files[new_index] = image_files->at(old_index);
    calibrated_intrinsics[new_index] = intrinsics->at(old_index);
    ++new_index;
  }

  *image_files = calibrated_image_files;
  *intrinsics = calibrated_intrinsics;

  for (auto& match : *matches) {
    const int new_image1_index =
        theia::FindOrDie(old_to_new_index, match.image1_index);
    const int new_image2_index =
        theia::FindOrDie(old_to_new_index, match.image2_index);
    match.image1_index = new_image1_index;
    match.image2_index = new_image2_index;
    match.twoview_info.focal_length_1 =
        intrinsics->at(new_image1_index).focal_length.value;
    match.twoview_info.focal_length_2 =
        intrinsics->at(new_image2_index).focal_length.value;
  }
}

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> image_files;
  std::vector<CameraIntrinsicsPrior> intrinsics;
  std::vector<ImagePairMatch> matches;

  // Read lists file.
  LOG(INFO) << "Readings lists file.";
  CHECK(ReadListsFile(FLAGS_list_file, &image_files));
  intrinsics.resize(image_files.size());

  // Read connected components.
  LOG(INFO) << "Reading connected components.";
  std::unordered_set<int> valid_images;
  CHECK(ReadCC(FLAGS_cc_file, &valid_images));

  // Read tracks.
  LOG(INFO) << "Reading tracks.";
  std::unordered_map<int, std::unordered_map<int, Eigen::Vector2d> > coords;
  std::unordered_map<int, std::vector<IntPair> > tracks;
  CHECK(ReadTracks(FLAGS_tracks_file, &coords, &tracks));

  // Read coords. Intrinsics are also provided with the coords file.
  LOG(INFO) << "Reading feature coordinates.";
  CHECK(ReadCoords(FLAGS_coords_file, valid_images, &intrinsics, &coords));

  // Read epipolar geometries.
  LOG(INFO) << "Reading epipolar geometries.";
  CHECK(ReadEGs(FLAGS_eg_file, &matches));

  LOG(INFO) << "Converting tracks to feature correspondences.";
  ConvertTracksToFeatureCorrespondences(coords, tracks, &matches);
  coords.clear();
  tracks.clear();

  LOG(INFO) << "Removing uncalibrated images and reindexing.";
  RemoveUncalibratedImageAndReindex(valid_images,
                                    &image_files,
                                    &intrinsics,
                                    &matches);

  LOG(INFO) << "Writing out matches.";
  CHECK(theia::WriteMatchesAndGeometry(FLAGS_output_match_file,
                                       image_files,
                                       intrinsics,
                                       matches))
      << "Could not write match file";

  return 0;
}
