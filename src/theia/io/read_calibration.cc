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

#include "theia/io/read_calibration.h"

#include <glog/logging.h>
#include <fstream>  // NOLINT
#include <iostream>  // NOLINT
#include <string>
#include <unordered_map>

#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {

bool ReadCalibration(const std::string& calibration_file,
                     std::unordered_map<std::string, CameraIntrinsicsPrior>*
                         camera_intrinsics_prior) {
  std::ifstream ifs(calibration_file.c_str(), std::ios::in);
  if (!ifs.is_open()) {
    LOG(ERROR) << "Cannot read the list file from " << calibration_file;
    return false;
  }

  while (!ifs.eof()) {
    // Read in the filename.
    std::string filename;
    ifs >> filename;
    if (filename.length() == 0) {
      break;
    }

    // Read camera_intrinsics_prior.
    CameraIntrinsicsPrior temp_camera_intrinsics_prior;
    temp_camera_intrinsics_prior.focal_length.is_set = true;
    ifs >> temp_camera_intrinsics_prior.focal_length.value;

    temp_camera_intrinsics_prior.principal_point[0].is_set = true;
    ifs >> temp_camera_intrinsics_prior.principal_point[0].value;
    temp_camera_intrinsics_prior.image_width =
        2.0 * temp_camera_intrinsics_prior.principal_point[0].value;

    temp_camera_intrinsics_prior.principal_point[1].is_set = true;
    ifs >> temp_camera_intrinsics_prior.principal_point[1].value;
    temp_camera_intrinsics_prior.image_height =
        2.0 * temp_camera_intrinsics_prior.principal_point[1].value;

    temp_camera_intrinsics_prior.aspect_ratio.is_set = true;
    ifs >> temp_camera_intrinsics_prior.aspect_ratio.value;

    temp_camera_intrinsics_prior.skew.is_set = true;
    ifs >> temp_camera_intrinsics_prior.skew.value;

    temp_camera_intrinsics_prior.radial_distortion[0].is_set = true;
    ifs >> temp_camera_intrinsics_prior.radial_distortion[0].value;

    temp_camera_intrinsics_prior.radial_distortion[1].is_set = true;
    ifs >> temp_camera_intrinsics_prior.radial_distortion[1].value;

    (*camera_intrinsics_prior)[filename] = temp_camera_intrinsics_prior;
  }
  return true;
}

}  // namespace theia
