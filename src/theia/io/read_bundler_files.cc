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

#include "theia/io/read_bundler_files.h"

#include <Eigen/Core>
#include <glog/logging.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>   // NOLINT
#include <iostream>  // NOLINT
#include <string>
#include <utility>
#include <vector>

#include "theia/io/bundler_file_reader.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/pinhole_camera_model.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"
#include "theia/util/filesystem.h"
#include "theia/util/map_util.h"

namespace theia {

namespace {

// // Description of the list files from the Big SfM website:
// // http://www.cs.cornell.edu/projects/p2f/README_Dubrovnik6K.txt
// //
// // A. List files (list.db.txt, list.query.txt).
// //      List files specify filenames to images in jpg format, one per
// //      line (keep in mind that the actual jpg files are not distributed
// //      unless requested).  In addition, if the focal length of the image
// //      has been estimated from Exif tags, then that is also included.
// //
// //      Images without known focal length information are specified with
// //      a line with a single field, the image name.  Example:
// //        query/10970812@N05_2553027508.jpg
// //
// //      Images with known focal length information are specified with a
// //      line with three fields: the image name, a zero, and the Exif
// //      focal length.  (The second field is always zero but may change in
// //      future datasets.)  Example:
// //        query/11289373@N03_2733280477.jpg 0 1280.00000
// //
// // NOTE: We set the exif focal length to zero if it is not available (since 0 is
// // never a valid focal length).
// bool ReadListsFile(const std::string& list_filename,
//                    Reconstruction* reconstruction) {
//   std::ifstream ifs(list_filename.c_str(), std::ios::in);
//   if (!ifs.is_open()) {
//     LOG(ERROR) << "Cannot read the list file from " << list_filename;
//     return false;
//   }

//   const char space = static_cast<char>(' ');
//   while (!ifs.eof()) {
//     // Read in the filename.
//     std::string filename, truncated_filename;
//     ifs >> filename;
//     if (filename.length() == 0) {
//       break;
//     }
//     CHECK(theia::GetFilenameFromFilepath(filename, true, &truncated_filename));
//     const ViewId view_id = reconstruction->AddView(truncated_filename);
//     CHECK_NE(view_id, kInvalidViewId)
//         << "View " << truncated_filename << " could not be added.";

//     // Check to see if the exif focal length is given.
//     double focal_length = 0;
//     if (ifs.peek() == space) {
//       int temp;
//       ifs >> temp;
//       ifs >> focal_length;
//     }

//     if (focal_length != 0) {
//       reconstruction->MutableView(view_id)
//           ->MutableCameraIntrinsicsPrior()
//           ->focal_length.value[0] = focal_length;
//       reconstruction->MutableView(view_id)
//           ->MutableCameraIntrinsicsPrior()
//           ->focal_length.is_set = true;
//     }
//   }
//   return true;
// }

bool AddViewsToReconstruction(const BundlerFileReader& reader,
                              Reconstruction* reconstruction) {
  if (reader.img_entries().empty()) {
    return false;
  }

  std::string truncated_filename;
  for (const ListImgEntry& entry : reader.img_entries()) {
    CHECK(theia::GetFilenameFromFilepath(
        entry.filename, true, &truncated_filename));
    const ViewId view_id = reconstruction->AddView(truncated_filename);
    CHECK_NE(view_id, kInvalidViewId)
        << "View " << truncated_filename << " could not be added.";
    // Set the focal length.
    if (entry.focal_length > 0.0f) {
      reconstruction->MutableView(view_id)
          ->MutableCameraIntrinsicsPrior()
          ->focal_length.value[0] = entry.focal_length;
      reconstruction->MutableView(view_id)
          ->MutableCameraIntrinsicsPrior()
          ->focal_length.is_set = true;
    }
  }

  return true;
}

bool AddCamerasToReconstruction(const BundlerFileReader& reader,
                                Reconstruction* reconstruction,
                                std::unordered_set<ViewId>* views_to_remove) {
  // Populate camera parameters.
  static const Eigen::Matrix3d bundler_to_theia =
      Eigen::Vector3d(1.0, -1.0, -1.0).asDiagonal();

  const int num_cameras = reader.NumCameras();
  CHECK_EQ(num_cameras, reconstruction->NumViews())
      << "The number of cameras in the lists file is not equal to the number "
      "of cameras in the bundle file. Data is corrupted!";

  // Read in the camera params.
  const std::vector<BundlerCamera>& bundler_cameras = reader.cameras();
  for (int i = 0; i < num_cameras; i++) {
    const BundlerCamera& bundler_camera = bundler_cameras[i];
    reconstruction->MutableView(i)->SetEstimated(true);
    Camera* camera = reconstruction->MutableView(i)->MutableCamera();
    camera->SetCameraIntrinsicsModelType(CameraIntrinsicsModelType::PINHOLE);

    // Do not consider this view if an invalid focal length is present.
    if (bundler_camera.focal_length <= 0.0) {
      views_to_remove->insert(i);
    } else {
      camera->SetFocalLength(bundler_camera.focal_length);
    }

    camera->MutableCameraIntrinsics()->SetParameter(
        PinholeCameraModel::RADIAL_DISTORTION_1, bundler_camera.radial_coeff_1);
    camera->MutableCameraIntrinsics()->SetParameter(
        PinholeCameraModel::RADIAL_DISTORTION_2, bundler_camera.radial_coeff_2);

    // These cameras (and the features below) already have the principal point
    // removed.
    camera->SetPrincipalPoint(0, 0);

    const Eigen::Matrix3d rotation =
        bundler_to_theia * bundler_camera.rotation;
    const Eigen::Vector3d translation = bundler_to_theia * bundler_camera.translation;

    const Eigen::Vector3d position = -rotation.transpose() * translation;
    camera->SetPosition(position);
    camera->SetOrientationFromRotationMatrix(rotation);

    if ((i + 1) % 100 == 0 || i == num_cameras - 1) {
      std::cout << "\r Loading parameters for camera " << i + 1 << " / "
                << num_cameras << std::flush;
    }
  }
  std::cout << std::endl;

  return true;
}

int AddTracksToReconstruction(
    const BundlerFileReader& reader,
    const  std::unordered_set<ViewId>& views_to_remove,
    Reconstruction* reconstruction) {
  // Read in each 3D point and correspondences.
  int num_invalid_tracks = 0;
  const int num_points = reader.NumPoints();
  const std::vector<BundlerPoint>& points = reader.points();
  for (int i = 0; i < num_points; i++) {
    const BundlerPoint point = points[i];
    const Eigen::Vector3d& position = point.position;
    const Eigen::Vector3d& color = point.color;
    const int num_views = point.view_list.size();

    // Reserve the view list for this 3D point.
    std::vector<std::pair<ViewId, Feature> > track;
    track.reserve(num_views);
    for (int j = 0; j < num_views; j++) {
      // TODO(vfragoso): Should we store SIFT indices? It is useful for
      // img-based localization.
      const FeatureInfo& feature_info = point.view_list[j];

      // NOTE: We flip the pixel directions to compensate for Bundlers different
      // coordinate system in images.
      const Feature feature(feature_info.kpt_x, -feature_info.kpt_y);

      // Push the sift key correspondence to the view list if the view is valid.
      if (!ContainsKey(views_to_remove, feature_info.camera_index)) {
        track.emplace_back(feature_info.camera_index, feature);
      }
    }

    // Do not add the track if it is underconstrained.
    if (track.size() < 2 || position.squaredNorm() == 0.0) {
      continue;
    }

    const TrackId track_id = reconstruction->AddTrack(track);
    if (track_id == kInvalidTrackId) {
      ++num_invalid_tracks;
      continue;
    }

    Track* mutable_track = reconstruction->MutableTrack(track_id);
    mutable_track->SetEstimated(true);
    *mutable_track->MutablePoint() = position.homogeneous();
    *mutable_track->MutableColor() = color.cast<uint8_t>();

    if ((i + 1) % 100 == 0 || i == num_points - 1) {
      std::cout << "\r Loading 3D points " << i + 1 << " / " << num_points
                << std::flush;
    }
  }

  return num_invalid_tracks;
}

}  // namespace

// The bundle files contain the estimated scene and camera geometry have the
// following format:
//     # Bundle file v0.3
//     <num_cameras> <num_points>   [two integers]
//     <camera1>
//     <camera2>
//        ...
//     <cameraN>
//     <point1>
//     <point2>
//        ...
//     <pointM>
// Each camera entry <cameraI> contains the estimated camera intrinsics and
// extrinsics, and has the form:
//     <f> <k1> <k2>   [the focal length, followed by two radial distortion
//                      coeffs]
//     <R>             [a 3x3 matrix representing the camera rotation]
//     <t>             [a 3-vector describing the camera translation]
// The cameras are specified in the order they appear in the list of images.
//
// Each point entry has the form:
//     <position>      [a 3-vector describing the 3D position of the point]
//     <color>         [a 3-vector describing the RGB color of the point]
//     <view list>     [a list of views the point is visible in]
//
// The view list begins with the length of the list (i.e., the number of cameras
// the point is visible in). The list is then given as a list of quadruplets
// <camera> <key> <x> <y>, where <camera> is a camera index, <key> the index of
// the SIFT keypoint where the point was detected in that camera, and <x> and
// <y> are the detected positions of that keypoint. Both indices are 0-based
// (e.g., if camera 0 appears in the list, this corresponds to the first camera
// in the scene file and the first image in "list.txt"). The pixel positions are
// floating point numbers in a coordinate system where the origin is the center
// of the image, the x-axis increases to the right, and the y-axis increases
// towards the top of the image. Thus, (-w/2, -h/2) is the lower-left corner of
// the image, and (w/2, h/2) is the top-right corner (where w and h are the
// width and height of the image).
bool ReadBundlerFiles(const std::string& lists_file,
                      const std::string& bundle_file,
                      Reconstruction* reconstruction) {
  CHECK_EQ(reconstruction->NumViews(), 0)
      << "An empty reconstruction must be provided to load a bundler dataset.";
  CHECK_EQ(reconstruction->NumTracks(), 0)
      << "An empty reconstruction must be provided to load a bundler dataset.";

  // Parse the bundler files.
  BundlerFileReader bundler_file_reader(lists_file, bundle_file);
  if (!bundler_file_reader.ParseListsFile()) {
    LOG(ERROR) << "Could not read the lists file from " << lists_file;
    return false;
  }

  if (!bundler_file_reader.ParseBundleFile()) {
    LOG(ERROR) << "Could not parse the bundler file from " << bundle_file;
    return false;
  }

  // Populate views in the reconstruction.
  CHECK(AddViewsToReconstruction(
      bundler_file_reader, CHECK_NOTNULL(reconstruction)));

  // Populate cameras to reconstruction.
  std::unordered_set<ViewId> views_to_remove;
  CHECK(AddCamerasToReconstruction(
      bundler_file_reader, reconstruction, &views_to_remove));

  // Populate tracks to reconstruction.
  const int num_invalid_tracks = AddTracksToReconstruction(
      bundler_file_reader, views_to_remove, reconstruction);

  // Remove any invalid views.
  for (const ViewId view_to_remove : views_to_remove) {
    reconstruction->RemoveView(view_to_remove);
  }

  const int num_cameras = bundler_file_reader.NumCameras();
  if (views_to_remove.size() > 0) {
    LOG(INFO) << "Could not add " << views_to_remove.size() << " out of "
              << num_cameras << " views due to invalid camera parameters.";
  }

  const int num_points = bundler_file_reader.NumPoints();
  if (num_invalid_tracks > 0) {
    LOG(INFO) << "Could not load " << num_invalid_tracks
              << " invalid tracks out of " << num_points
              << " total tracks from Bundler file. This typically occurs when "
                 "tracks contain multiple observations from the same image.";
  }

  return true;
}

}  // namespace theia
