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
// Author: Victor Fragoso (victor.fragoso@mail.wvu.edu)

#include "theia/io/bundler_file_reader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <glog/logging.h>

#include "theia/sfm/reconstruction.h"

namespace theia {
namespace {
bool ReadHeader(FILE* in, int* num_cameras, int* num_points) {
  // Read the comment.
  char* line = nullptr;
  size_t length = 0;
  if ((getline(&line, &length, in) == -1) || line[0] != '#') {
    free(line);
    return false;
  }
  VLOG(3) << "Comment: " << line;
  free(line);
  // Read number of points and cameras.
  if (fscanf(in, "%d %d",
             CHECK_NOTNULL(num_cameras),
             CHECK_NOTNULL(num_points)) != 2 ||
      *num_cameras < 1 || *num_points < 1) {
    return false;
  }
  VLOG(3) << "Num cameras to read: " << *num_cameras;
  VLOG(3) << "Num points to read: " << *num_points;
  return true;
}

bool ReadCamera(FILE* in, BundlerCamera* camera) {
  // Read focal length and radial distortion coeffs.
  if (fscanf(in, "%f %f %f",
             &camera->focal_length,
             &camera->radial_coeff_1,
             &camera->radial_coeff_2) != 3) {
    VLOG(3) << "Unable to read focal length and radial distortion coeffs.";
    return false;
  }
  // Read rotation matrix.
  Eigen::Matrix3d& rotation = camera->rotation;
  if (fscanf(in, "%lf %lf %lf",
             &rotation(0, 0), &rotation(0, 1), &rotation(0, 2)) != 3) {
    VLOG(3) << "Unable to read first row of rotation matrix.";
    return false;
  }
  if (fscanf(in, "%lf %lf %lf",
             &rotation(1, 0), &rotation(1, 1), &rotation(1, 2)) != 3) {
    VLOG(3) << "Unable to read second row of rotation matrix.";
    return false;
  }
  if (fscanf(in, "%lf %lf %lf",
             &rotation(2, 0), &rotation(2, 1), &rotation(2, 2)) != 3) {
    VLOG(3) << "Unable to read third row of rotation matrix.";
    return false;
  }
  // Read position.
  Eigen::Vector3d& translation = camera->translation;
  if (fscanf(in, "%lf %lf %lf",
             &translation(0), &translation(1), &translation(2)) != 3) {
    VLOG(3) << "Unable to read camera translation.";
    return false;
  }
  return true;
}

bool ReadCameras(const int num_cameras,
                 FILE* in,
                 std::vector<BundlerCamera>* cameras) {
  CHECK_NOTNULL(cameras)->reserve(num_cameras);
  for (int i = 0; i < num_cameras; ++i) {
    cameras->emplace_back();
    if (!ReadCamera(in, &cameras->back())) {
      return false;
    }
  }
  VLOG(3) << "First camera in file: ";
  VLOG(3) << "Focal length: " << cameras->front().focal_length;
  VLOG(3) << "Dist coeff 1: " << cameras->front().radial_coeff_1;
  VLOG(3) << "Dist coeff 2: " << cameras->front().radial_coeff_2;
  VLOG(3) << "Rotation: \n" << cameras->front().rotation;
  VLOG(3) << "Translation: \n" << cameras->front().translation.transpose();
  VLOG(3) << "Last camera in file: ";
  VLOG(3) << "Focal length: " << cameras->back().focal_length;
  VLOG(3) << "Dist coeff 1: " << cameras->back().radial_coeff_1;
  VLOG(3) << "Dist coeff 2: " << cameras->back().radial_coeff_2;
  VLOG(3) << "Rotation: \n" << cameras->back().rotation;
  VLOG(3) << "Translation: \n" << cameras->back().translation.transpose();
  return true;
}

bool ReadViewList(FILE* in, std::vector<FeatureInfo>* view_list) {
  // Read number of views.
  int num_views = 0;
  if (fscanf(in, "%d", &num_views) != 1) {
    VLOG(3) << "Unable to read number of views for point.";
  }
  VLOG(3) << "Num. views to read: " << num_views;
  view_list->resize(num_views);
  float entry = 0;
  for (int i = 0; i < view_list->size(); ++i) {
    // Camera index.
    if (fscanf(in, "%f", &entry) != 1) {
      return false;
    }
    (*view_list)[i].camera_index = static_cast<int>(entry);
    // Sift index.
    if (fscanf(in, "%f", &entry) != 1) {
      return false;
    }
    (*view_list)[i].sift_index = static_cast<int>(entry);
    // Kpt x.
    if (fscanf(in, "%f", &entry) != 1) {
      return false;
    }
    (*view_list)[i].kpt_x = static_cast<int>(entry);
    // Kpt y.
    if (fscanf(in, "%f", &entry) != 1) {
      return false;
    }
    (*view_list)[i].kpt_y = static_cast<int>(entry);
  }
  return true;
}

bool ReadPoint(FILE* in, BundlerPoint* point) {
  // Read position.
  Eigen::Vector3d& position = point->position;
  if (fscanf(in, "%lf %lf %lf",
             &position(0), &position(1), &position(2)) != 3) {
      VLOG(3) << "Unable to read point position. ";
    return false;
  }
  // Read color.
  Eigen::Vector3d& color = point->color;
  if (fscanf(in, "%lf %lf %lf",
             &color(0), &color(1), &color(2)) != 3) {
    VLOG(3) << "Unable to read point color.";
    return false;
  }
  // Read the view list.
  if (!ReadViewList(in, &point->view_list)) {
    VLOG(3) << "Unable to read view list.";
    return false;
  }
  return true;
}

bool ReadPoints(const int num_points,
                FILE* in,
                std::vector<BundlerPoint>* points) {
  CHECK_NOTNULL(points)->reserve(num_points);
  for (int i = 0; i < num_points; ++i) {
    points->emplace_back();
    if (!ReadPoint(in, &points->back())) {
      return false;
    }
  }
  return true;
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
bool BundlerFileReader::ParseBundleFile() {
  FILE* in = fopen(bundler_filepath_.c_str(), "r");
  if (in == nullptr) {
    LOG(INFO) << "Could not open: " << bundler_filepath_;
    return false;
  }
  // Read Header.
  int num_cameras;
  int num_points;
  if (!ReadHeader(in, &num_cameras, &num_points)) {
    fclose(in);
    VLOG(3) << "Unable to read header.";
    return false;
  }
  // Read Cameras.
  if (!ReadCameras(num_cameras, in, &cameras_)) {
    fclose(in);
    VLOG(3) << "Unable to read cameras.";
    return false;
  }
  // Read Points.
  if (!ReadPoints(num_points, in, &points_)) {
    fclose(in);
    VLOG(3) << "Unable to read points.";
    return false;
  }
  // Close the file.
  fclose(in);
  bundler_file_parsed_ = true;
  return true;
}


// Description of the list files from the Big SfM website:
// http://www.cs.cornell.edu/projects/p2f/README_Dubrovnik6K.txt
//
// A. List files (list.db.txt, list.query.txt).
//      List files specify filenames to images in jpg format, one per
//      line (keep in mind that the actual jpg files are not distributed
//      unless requested).  In addition, if the focal length of the image
//      has been estimated from Exif tags, then that is also included.
//
//      Images without known focal length information are specified with
//      a line with a single field, the image name.  Example:
//        query/10970812@N05_2553027508.jpg
//
//      Images with known focal length information are specified with a
//      line with three fields: the image name, a zero, and the Exif
//      focal length.  (The second field is always zero but may change in
//      future datasets.)  Example:
//        query/11289373@N03_2733280477.jpg 0 1280.00000
//
// NOTE: We set the exif focal length to zero if it is not available (since 0 is
// never a valid focal length).
bool BundlerFileReader::ParseListsFile() {
  FILE* in = fopen(lists_filepath_.c_str(), "r");
  if (in == nullptr) {
    LOG(INFO) << "Could not open: " << bundler_filepath_;
    return false;
  }
  // Read line by line
  char* line = nullptr;
  size_t length = 0;
  img_entries_.reserve(1024);
  char filename[1024];
  while (getline(&line, &length, in) != -1) {
    img_entries_.emplace_back();
    ListImgEntry& entry = img_entries_.back();
    const int read_items = sscanf(line, "%s %f %f",
           filename, &entry.second_entry, &entry.focal_length);
    if (read_items != 1 && read_items !=3 ) {
      VLOG(3) << "Invalid line: " << line;
      return false;
    }
    entry.filename = filename;
  }
  free(line);
  fclose(in);
  lists_file_parsed_ = true;
  return true;
}

}  // namespace theia
