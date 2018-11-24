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

#ifndef THEIA_IO_BUNDLER_FILE_READER_H_
#define THEIA_IO_BUNDLER_FILE_READER_H_

#include <string>
#include <vector>

#include <Eigen/Core>

#include "theia/sfm/reconstruction.h"

namespace theia {

// Camera extrinsics and intrinsics.
struct BundlerCamera {
  // Rigid transformation from world to camera.
  Eigen::Vector3d translation;
  Eigen::Matrix3d rotation;
  // Estimated focal lengths and radial dist coeffs.
  float focal_length;
  float radial_coeff_1;
  float radial_coeff_2;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Feature2d information.
struct FeatureInfo {
  // Index of the camera in the lists file.
  int camera_index;
  // Index of the feature in the SIFT file.
  int sift_index;
  // Pixel position.
  int kpt_x;
  int kpt_y;
};

// Reconstructed 3D point.
struct BundlerPoint {
  // 3D point's position.
  Eigen::Vector3d position;
  // 3D point's color.
  Eigen::Vector3d color;
  // A list of where this poins is viewed.
  std::vector<FeatureInfo> view_list;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Entry in the lists.txt file.
struct ListImgEntry {
  // Image filename.
  std::string filename;
  // Second entry in the line.
  float second_entry;
  // Focal length.
  float focal_length;
};

class BundlerFileReader {
 public:
  // Params:
  //   lists_filepath  The filepath of the lists.txt file with the list of
  //     images.
  //   bundler_filepath  The filepath of the bundler output file.
  BundlerFileReader(const std::string& lists_filepath,
                    const std::string& bundler_filepath) :
      lists_filepath_(lists_filepath),
      bundler_filepath_(bundler_filepath),
      bundler_file_parsed_(false),
      lists_file_parsed_(false) {}
  virtual ~BundlerFileReader() {}

  // Parses the bundler output file: bundler_filepath_. Returns true upon
  // success, and false otherwise.
  bool ParseBundleFile();

  // Parses the lists.txt file. Returns true upon success, and false otherwise.
  bool ParseListsFile();

  // Returns the number of cameras contained in the bundler output file.
  int NumCameras() { return cameras_.size();}
  int NumCameras() const { return cameras_.size();}

  // Returns the number of reconstructed points in the bundler output file.
  int NumPoints() { return points_.size();}
  int NumPoints() const { return points_.size();}

  // Returns the number of entries in the lists.txt file.
  int NumListEntries() { return img_entries_.size();}
  int NumListEntries() const { return img_entries_.size();}

  // Getters and setters.
  const std::vector<BundlerCamera>& cameras() const {
    return cameras_;
  }

  const std::vector<BundlerPoint>& points() const {
    return points_;
  }

  std::vector<BundlerPoint>* mutable_points()  {
    return &points_;
  }

  const std::vector<ListImgEntry>& img_entries() const {
    return img_entries_;
  }
  
 private:
  // Lists filepath.
  const std::string lists_filepath_;
  // Bundler filepath.
  const std::string bundler_filepath_;
  // Bundler cameras.
  std::vector<BundlerCamera> cameras_;
  // Bundler points.
  std::vector<BundlerPoint> points_;
  // List image entries.
  std::vector<ListImgEntry> img_entries_;
  // Flags indicating that a file was parsed.
  bool bundler_file_parsed_;
  bool lists_file_parsed_;
};

}  // namespace theia

#endif  // THEIA_IO_BUNDLER_FILE_READER_H_
