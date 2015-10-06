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

#include "theia/io/sift_text_file.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>  // NOLINT
#include <iostream>  // NOLINT
#include <string>
#include <vector>

namespace theia {

// The sift key file has the following format:
//
// number_of_keypoints sift_descriptor_dimensions (both as ints)
// for each descriptor:
//   row col scale orientation (all as floats)
//   128 ints describing sift descriptor. Normalizing this 128-vector to unit
//     length will yield the true sift descriptor.
// NOTE: We use getline and strtof which should be much faster than letting the
// stream parse the string with operator >>.
bool ReadSiftKeyTextFile(const std::string& sift_key_file,
                         std::vector<Eigen::VectorXf>* descriptor,
                         std::vector<Keypoint>* keypoint) {
  CHECK_NOTNULL(descriptor)->clear();
  CHECK_NOTNULL(keypoint)->clear();

  FILE* fp = fopen(sift_key_file.c_str(), "r");
  int num_descriptors, len;

  if (fscanf(fp, "%d %d", &num_descriptors, &len) != 2) {
    printf("Invalid keypoint file\n");
    return 0;
  }

  CHECK_EQ(len, 128);

  descriptor->reserve(num_descriptors);
  keypoint->reserve(num_descriptors);

  for (int i = 0; i < num_descriptors; i++) {
    float x, y, scale, ori;

    if (fscanf(fp, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) {
      printf("Invalid keypoint file format.");
      return 0;
    }

    Keypoint kp(x, y, Keypoint::SIFT);
    kp.set_scale(scale);
    kp.set_orientation(ori);
    keypoint->push_back(kp);

    char buf[1024];
    Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> int_descriptor(128);
    uint8_t* p = int_descriptor.data();
    for (int line = 0; line < 6; line++) {
      CHECK_NOTNULL(fgets(buf, 1024, fp));
      sscanf(buf, "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu "
             "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
             p + 0, p + 1, p + 2, p + 3, p + 4, p + 5, p + 6, p + 7, p + 8,
             p + 9, p + 10, p + 11, p + 12, p + 13, p + 14, p + 15, p + 16,
             p + 17, p + 18, p + 19);

      p += 20;
    }
    CHECK_NOTNULL(fgets(buf, 1024, fp));
    sscanf(buf, "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu", p + 0, p + 1,
           p + 2, p + 3, p + 4, p + 5, p + 6, p + 7);

    Eigen::VectorXf float_descriptor = int_descriptor.cast<float>();
    float_descriptor /= 255.0;
    descriptor->push_back(float_descriptor);
  }

  fclose(fp);
  return true;
}

}  // namespace theia
