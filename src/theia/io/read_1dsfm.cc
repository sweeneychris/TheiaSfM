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

#include "theia/io/read_1dsfm.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>

#include <algorithm>
#include <fstream>  // NOLINT
#include <string>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "theia/sfm/find_common_tracks_in_views.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/sfm/view.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/util/filesystem.h"
#include "theia/util/map_util.h"

namespace theia {

class Input1DSFM {
 public:
  Input1DSFM(const std::string& dataset_directory,
             Reconstruction* reconstruction,
             ViewGraph* view_graph)
      : dataset_directory_(dataset_directory),
        reconstruction_(reconstruction),
        view_graph_(view_graph) {}

  bool ReadCC(std::unordered_set<int>* valid_image_index);
  bool ReadListsFile(const std::unordered_set<int>& valid_image_index);
  bool ReadTracks();
  bool ReadCoords();
  bool ReadEGs();

 private:
  bool ReadCoordsHeaderLine(const std::string& line,
                            ViewId* image_index,
                            int* num_keys);

  const std::string& dataset_directory_;
  Reconstruction* reconstruction_;
  ViewGraph* view_graph_;

  // Maps image id -> [ feature id -> feature coordinate].
  std::unordered_map<ViewId, std::vector<Feature> > feature_coordinates_;
};

// Reads the connected components file.
bool Input1DSFM::ReadCC(std::unordered_set<int>* valid_image_index) {
  const std::string cc_filename = dataset_directory_ + "/cc.txt";

  std::ifstream ifs(cc_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the cc file from " << cc_filename;
    return false;
  }

  while (!ifs.eof()) {
    int img_index;
    ifs >> img_index;
    valid_image_index->insert(img_index);
  }

  return true;
}

// Reads the list file and ignores the focal length since it can be recovered
// from the coords file.
bool Input1DSFM::ReadListsFile(
    const std::unordered_set<int>& valid_image_index) {
  const std::string list_filename = dataset_directory_ + "/list.txt";
  std::ifstream ifs(list_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the list file from " << list_filename;
    return false;
  }

  const char space = static_cast<char>(' ');
  while (!ifs.eof()) {
    // Read in the filename.
    std::string filename, truncated_filename;

    ifs >> filename;
    if (filename.length() == 0) {
      break;
    }
    CHECK(theia::GetFilenameFromFilepath(filename, false, &truncated_filename));
    const ViewId view_id = reconstruction_->AddView(truncated_filename);
    CHECK_NE(view_id, kInvalidViewId);

    // Check to see if the exif focal length is given.
    double focal_length = 0;
    if (ifs.peek() == space) {
      int temp;
      ifs >> temp;
      ifs >> focal_length;
    }

    // If the view is not in the connected component remove it. Adding it first
    // then removing it allows for our ViewIds to stay in sync with the index of
    // the image in hte lists file. This will simplify the indexing throughout
    // the entire program.
    if (!ContainsKey(valid_image_index, view_id)) {
      reconstruction_->RemoveView(view_id);
      continue;
    }

    if (focal_length != 0) {
      reconstruction_->MutableView(view_id)
          ->MutableCameraIntrinsicsPrior()
          ->focal_length.value = focal_length;
      reconstruction_->MutableView(view_id)
          ->MutableCameraIntrinsicsPrior()
          ->focal_length.is_set = true;
      LOG(INFO) << "Adding image " << truncated_filename
                << " with focal length: " << focal_length;
    } else {
      LOG(INFO) << "Adding image " << truncated_filename
                << " with focal length: UNKNOWN";
    }
  }
  return true;
}

bool Input1DSFM::ReadCoordsHeaderLine(const std::string& line,
                                      ViewId* view_id,
                                      int* num_keys) {
  float principal_point_x, principal_point_y, focal_length;
  char name[256];
  sscanf(line.c_str(),
         "#index = %d, name = %s keys = %d, px = %f, py = %f, focal = %f",
         view_id,
         name,
         num_keys,
         &principal_point_x,
         &principal_point_y,
         &focal_length);

  View* view = reconstruction_->MutableView(*view_id);
  if (view == nullptr) {
    return false;
  }

  // Set the metadata.
  CameraIntrinsicsPrior* prior = view->MutableCameraIntrinsicsPrior();
  prior->image_width = principal_point_x * 2.0;
  prior->image_width = principal_point_y * 2.0;
  prior->principal_point[0].value = principal_point_x;
  prior->principal_point[0].is_set = true;
  prior->principal_point[1].value = principal_point_y;
  prior->principal_point[1].is_set = true;
  return true;
}

// Reads the coords file. Only the coords with a valid track in the connected
// component are kept.
bool Input1DSFM::ReadCoords() {
  const std::string coords_filename = dataset_directory_ + "/coords.txt";
  std::ifstream ifs(coords_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the coords file from " << coords_filename;
    return false;
  }

  feature_coordinates_.reserve(reconstruction_->NumViews());
  while (!ifs.eof()) {
    std::string line;
    std::getline(ifs, line);
    if (ifs.eof()) {
      break;
    }

    int num_keys;
    ViewId view_id;
    // If the image is not in the connected component then do not read it.
    if (!ReadCoordsHeaderLine(line, &view_id, &num_keys)) {
      for (int i = 0; i < num_keys; i++) {
        std::getline(ifs, line);
      }
      continue;
    }

    auto& features = feature_coordinates_[view_id];
    features.reserve(num_keys);

    Eigen::Vector2d keypoint;
    for (int i = 0; i < num_keys; i++) {
      std::getline(ifs, line);
      sscanf(line.c_str(), "%*d %lf %lf", &keypoint[0], &keypoint[1]);
      features.emplace_back(keypoint);
    }
  }

  return true;
}

bool Input1DSFM::ReadTracks() {
  const std::string tracks_filename = dataset_directory_ + "/tracks.txt";
  std::ifstream ifs(tracks_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the coords file from " << tracks_filename;
    return false;
  }

  // Read number of tracks.
  int num_tracks;
  ifs >> num_tracks;

  for (int i = 0; i < num_tracks; i++) {
    int num_features;
    ifs >> num_features;

    std::vector<std::pair<ViewId, Feature> > track;
    track.reserve(num_features);
    int feature_id;
    ViewId view_id;
    for (int j = 0; j < num_features; j++) {
      ifs >> view_id;
      ifs >> feature_id;

      const auto& features = FindOrDie(feature_coordinates_, view_id);
      const Feature& feature = features[feature_id];
      track.emplace_back(view_id, feature);
    }
    CHECK_NE(reconstruction_->AddTrack(track), kInvalidTrackId);
  }

  return true;
}

// Reads the epipolar geometry files.
bool Input1DSFM::ReadEGs() {
  const std::string eg_filename = dataset_directory_ + "/EGs.txt";
  std::ifstream ifs(eg_filename.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the EG file from " << eg_filename;
    return false;
  }

  const Eigen::Matrix3d bundler_to_theia =
      Eigen::Vector3d(1.0, -1.0, -1.0).asDiagonal();
  while (!ifs.eof()) {
    TwoViewInfo info;
    ViewId view_id1, view_id2;
    ifs >> view_id1;
    ifs >> view_id2;

    // The rotation defines the camera 2 to camera 1 transformation in row-major
    // order). We want a camera 1 to camera 2 transformation so we read in the
    // transpose (i.e., column-major order).
    Eigen::Matrix3d rotation;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        ifs >> rotation(i, j);
      }
    }

    rotation = bundler_to_theia * rotation.transpose() * bundler_to_theia;

    // Convert to angle axis.
    ceres::RotationMatrixToAngleAxis(rotation.data(), info.rotation_2.data());

    // Assign the position.
    ifs >> info.position_2[0];
    ifs >> info.position_2[1];
    ifs >> info.position_2[2];

    // Add the focal lengths. If they are known from EXIF, add that value
    // otherwise add a focal length guess correspdonding to a median viewing
    // angle.
    const CameraIntrinsicsPrior prior1 =
        reconstruction_->View(view_id1)->CameraIntrinsicsPrior();
    const CameraIntrinsicsPrior prior2 =
        reconstruction_->View(view_id2)->CameraIntrinsicsPrior();
    if (prior1.focal_length.is_set) {
      info.focal_length_1 = prior1.focal_length.value;
    } else {
      info.focal_length_1 = 1.2 * prior1.principal_point[0].value;
    }

    if (prior2.focal_length.is_set) {
      info.focal_length_2 = prior2.focal_length.value;
    } else {
      info.focal_length_2 = 1.2 * prior2.principal_point[0].value;
    }

    // Add the number of inliers.
    const std::vector<ViewId> views = {view_id1, view_id2};
    const std::vector<TrackId> common_tracks =
        FindCommonTracksInViews(*reconstruction_, views);
    info.num_verified_matches = common_tracks.size();

    // Add the match to the output.
    if (reconstruction_->View(view_id1) != nullptr &&
        reconstruction_->View(view_id2) != nullptr) {
      view_graph_->AddEdge(view_id1, view_id2, info);
    }
  }
  return true;
}

bool Read1DSFM(const std::string& dataset_directory,
               Reconstruction* reconstruction,
               ViewGraph* view_graph) {
  CHECK_NOTNULL(reconstruction);
  CHECK_NOTNULL(view_graph);

  Input1DSFM input_reader(dataset_directory, reconstruction, view_graph);

  LOG(INFO) << "Reading connected components.";
  std::unordered_set<int> valid_images;
  if (!input_reader.ReadCC(&valid_images)) {
    return false;
  }

  LOG(INFO) << "Readings lists file.";
  if (!input_reader.ReadListsFile(valid_images)) {
    return false;
  }

  LOG(INFO) << "Reading feature coordinates.";
  if (!input_reader.ReadCoords()) {
    return false;
  }

  LOG(INFO) << "Reading tracks.";
  if (!input_reader.ReadTracks()) {
    return false;
  }

  LOG(INFO) << "Reading epipolar geometries.";
  if (!input_reader.ReadEGs()) {
    return false;
  }

  return true;
}

}  // namespace theia
