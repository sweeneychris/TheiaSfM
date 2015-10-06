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

#include "theia/io/sift_binary_file.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <stdint.h>

#include <cstdlib>
#include <fstream>  // NOLINT
#include <iostream>  // NOLINT
#include <string>
#include <vector>

#include "theia/image/keypoint_detector/keypoint.h"

namespace theia {

// The sift key file has the following format:
//
// number_of_keypoints sift_descriptor_dimensions (both as ints)
// for each descriptor:
//   row col scale orientation (all as floats)
//   128 ints describing sift descriptor. Normalizing this 128-vector to unit
//     length will yield the float sift descriptor.
bool ReadSiftKeyBinaryFile(const std::string& sift_key_file,
                           std::vector<Eigen::VectorXf>* descriptor,
                           std::vector<Keypoint>* keypoint) {
  CHECK_NOTNULL(descriptor)->clear();
  CHECK_NOTNULL(keypoint)->clear();

  std::ifstream ifs(sift_key_file.c_str(), std::ios::in | std::ios::binary);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Could not read the sift key binary file from "
               << sift_key_file;
    return false;
  }

  // Number of descriptors and the length of the descriptors.
  int num_descriptors, len;
  ifs.read(reinterpret_cast<char*>(&num_descriptors), sizeof(num_descriptors));
  ifs.read(reinterpret_cast<char*>(&len), sizeof(len));
  CHECK_EQ(len, 128);

  descriptor->reserve(num_descriptors);
  keypoint->reserve(num_descriptors);

  for (int i = 0; i < num_descriptors; i++) {
    // Keypoint params = y, x, scale, orientation.
    float keypoint_params[4];
    ifs.read(reinterpret_cast<char*>(keypoint_params),
             4 * sizeof(keypoint_params[0]));

    Keypoint kp(keypoint_params[1], keypoint_params[0], Keypoint::SIFT);
    kp.set_scale(keypoint_params[2]);
    kp.set_orientation(keypoint_params[3]);
    keypoint->push_back(kp);

    Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> int_desc(128);
    ifs.read(reinterpret_cast<char*>(int_desc.data()),
             int_desc.size() * sizeof(int_desc[0]));
    Eigen::VectorXf float_descriptor = int_desc.cast<float>();
    float_descriptor /= 255.0;
    descriptor->push_back(float_descriptor);
  }
  ifs.close();

  return true;
}

// Outputs the SIFT features in the same format as Lowe's sift key files, but
// stores it as a binary file for faster loading.
bool WriteSiftKeyBinaryFile(
    const std::string& output_sift_key_file,
    const std::vector<Eigen::VectorXf>& descriptor,
    const std::vector<Keypoint>& keypoint) {
  CHECK_EQ(descriptor.size(), keypoint.size());

  std::ofstream ofs(output_sift_key_file.c_str(),
                    std::ios::out | std::ios::binary);
  if (!ofs.is_open()) {
    LOG(ERROR) << "Could not write the sift key binary file to "
               << output_sift_key_file;
    return false;
  }

  // Output number of descriptors and descriptor length.
  const int num_descriptors = descriptor.size();
  const int len = 128;
  ofs.write(reinterpret_cast<const char*>(&num_descriptors),
            sizeof(num_descriptors));
  ofs.write(reinterpret_cast<const char*>(&len), sizeof(len));

  // Output the sift descriptors.
  for (int i = 0; i < num_descriptors; i++) {
    // Output the keypoint information (x, y, scale, orientation).
    float kp[4];
    kp[0] = keypoint[i].y();
    kp[1] = keypoint[i].x();
    kp[2] = keypoint[i].scale();
    kp[3] = keypoint[i].orientation();
    ofs.write(reinterpret_cast<const char*>(kp), 4 * sizeof(kp[0]));

    // Output the descriptor.
    const Eigen::VectorXf float_scaled_desc = descriptor[i] * 255.0;
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> int_desc =
        float_scaled_desc.cast<uint8_t>();
    ofs.write(reinterpret_cast<const char*>(int_desc.data()),
              int_desc.size() * sizeof(int_desc[0]));
  }
  ofs.close();
  return true;
}

}  // namespace theia
