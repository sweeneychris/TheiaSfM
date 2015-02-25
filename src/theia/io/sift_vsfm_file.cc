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
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu), Benjamin Nuernberger
// (bnuernberger@cs.ucsb.edu).


#include "theia/io/sift_vsfm_file.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <cstdlib>
#include <fstream>  // NOLINT
#include <iostream>  // NOLINT
#include <string>
#include <vector>

#include "theia/image/keypoint_detector/keypoint.h"

namespace theia {

// These functions are modeled after the code found in the Viewpoint Invariant
// Patch code at http://ccwu.me/code.html

// The VisualSFM sift file has the following format
// (per http://ccwu.me/vsfm/doc.html#customize):
// 
// [Header] = int[5] = {name, version, npoint, 5, 128};
// name = ('S'+ ('I'<<8)+('F'<<16)+('T'<<24));
// version = ('V'+('4'<<8)+('.'<<16)+('0'<<24)); 
//        or ('V'+('5'<<8)+('.'<<16)+('0'<<24)) if containing color info
// npoint = number of features.
// 
// [Location Data] is a npoint x 5 float matrix and each row is:
//        [x, y, color, scale, orientation] 
// Color is read by casting the float to unsigned char[4].
// 
// [Descriptor Data] is a npoint x 128 unsigned char matrix.
// Note the feature descriptors are normalized to 512.
// 
// [EOF]  int eof_marker = (0xff+('E'<<8)+('O'<<16)+('F'<<24));
//
// Notes:
// * VisualSFM sorts the features in the order of decreasing scales.
//
bool ReadSiftVisualSFMFile(const std::string& input_sift_vsfm_file,
                           std::vector<Eigen::VectorXf>* descriptors,
                           std::vector<Keypoint>* keypoints) {
  CHECK_NOTNULL(descriptors)->clear();
  CHECK_NOTNULL(keypoints)->clear();

  std::ifstream ifs(input_sift_vsfm_file.c_str(), 
                    std::ios::in | std::ios::binary);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Could not read the VisualSFM sift key binary file from "
               << input_sift_vsfm_file;
    return false;
  }

  // Read header
  int name, version, num_descriptors, location_length, descriptor_length;
  ifs.read(reinterpret_cast<char*>(&name), sizeof(name));
  ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
  ifs.read(reinterpret_cast<char*>(&num_descriptors), sizeof(num_descriptors));
  ifs.read(reinterpret_cast<char*>(&location_length), sizeof(location_length));
  ifs.read(reinterpret_cast<char*>(&descriptor_length),
           sizeof(descriptor_length));

  CHECK_EQ(location_length, 5);
  CHECK_EQ(descriptor_length, 128);

  descriptors->reserve(num_descriptors);
  keypoints->reserve(num_descriptors);

  for (int i = 0; i < num_descriptors; i++) {
    // Keypoint params = y, x, color, scale, orientation.
    float keypoint_params[5];
    ifs.read(reinterpret_cast<char*>(keypoint_params),
             5 * sizeof(keypoint_params[0]));

    Keypoint kp(keypoint_params[1], keypoint_params[0], Keypoint::SIFT);
    // TODO: color?
    kp.set_scale(keypoint_params[3]);
    kp.set_orientation(keypoint_params[4]);
    keypoints->push_back(kp);
  }

  for (int i = 0; i < num_descriptors; i++) {
    Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> int_desc(128);
    ifs.read(reinterpret_cast<char*>(int_desc.data()),
             int_desc.size() * sizeof(int_desc[0]));
    Eigen::VectorXf float_descriptor = int_desc.cast<float>();
    float_descriptor /= 255.0;
    descriptors->push_back(float_descriptor);
  }
  ifs.close();

  return true;
}

// Outputs the SIFT features in the same format as VisualSFM's sift key files.
bool WriteSiftVisualSFMFile(
    const std::string& output_sift_vsfm_file,
    const std::vector<Eigen::VectorXf>& descriptors,
    const std::vector<Keypoint>& keypoints) {
  CHECK_EQ(descriptors.size(), keypoints.size());

  std::ofstream ofs(output_sift_vsfm_file.c_str(),
                    std::ios::out | std::ios::binary);
  if (!ofs.is_open()) {
    LOG(ERROR) << "Could not write the VisualSFM sift key binary file to "
               << output_sift_vsfm_file;
    return false;
  }

  // Output the header
  const int SIFT_NAME = ('S' + ('I'<<8) + ('F'<<16) + ('T'<<24));
  const int SIFT_VERSION_4 = ('V' + ('4'<<8) + ('.'<<16) + ('0'<<24));
  const int num_descriptors = descriptors.size();
  const int location_length = 5;
  const int descriptor_length = 128;
  ofs.write(reinterpret_cast<const char*>(&SIFT_NAME), sizeof(SIFT_NAME));
  ofs.write(reinterpret_cast<const char*>(&SIFT_VERSION_4),
            sizeof(SIFT_VERSION_4));
  ofs.write(reinterpret_cast<const char*>(&num_descriptors),
            sizeof(num_descriptors));
  ofs.write(reinterpret_cast<const char*>(&location_length),
            sizeof(location_length));
  ofs.write(reinterpret_cast<const char*>(&descriptor_length),
            sizeof(descriptor_length));

  // Output the keypoint information (x, y, scale, orientation).
  for (int i = 0; i < num_descriptors; i++) {
    float kp[5];
    kp[0] = keypoints[i].y();
    kp[1] = keypoints[i].x();
    kp[2] = 0; // TODO: color?
    kp[3] = keypoints[i].scale();
    kp[4] = keypoints[i].orientation();
    ofs.write(reinterpret_cast<const char*>(kp), 5 * sizeof(kp[0]));
  }

  // Output the sift descriptors.
  for (int i = 0; i < num_descriptors; i++) {
    const Eigen::VectorXf float_scaled_desc = descriptors[i] * 255.0;
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> int_desc =
        float_scaled_desc.cast<uint8_t>();
    ofs.write(reinterpret_cast<const char*>(int_desc.data()),
              int_desc.size() * sizeof(int_desc[0]));
  }
  ofs.close();
  return true;
}

}  // namespace theia
