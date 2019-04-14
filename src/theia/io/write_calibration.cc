// Copyright (C) 2019 The Regents of the University of California (Regents).
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

#include "theia/io/write_calibration.h"

#include <fstream>  // NOLINT
#include <string>
#include <unordered_map>

#include <glog/logging.h>

#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {
namespace {

bool WritePrior(const std::string& img_name,
                const CameraIntrinsicsPrior& prior,
                std::ofstream* out_ptr) {
  std::ofstream& out = *out_ptr;
  out << "{\"CameraIntrinsicsPrior\" : {\n";
  // Writing the attributes that must be set by default.
  out << "\"image_name\" : \"" << img_name << "\",\n";
  out << "\"width\" : " << prior.image_width << ",\n";
  out << "\"height\" : " << prior.image_height << ",\n";
  out << "\"camera_intrinsics_type\" : \""
      << prior.camera_intrinsics_model_type << "\"";
  // Place coma before any of the additional attributes if one is set.
  if (prior.focal_length.is_set) {
    out << ",\n";
    out << "\"focal_length\" : " << prior.focal_length.value[0];
  }
  if (prior.principal_point.is_set) {
    out << ",\n";
    out << "\"principal_point\" : ["
        << prior.principal_point.value[0] << ", "
        << prior.principal_point.value[1] << "]";
  }
  if (prior.aspect_ratio.is_set) {
    out << ",\n";
    out << "\"aspect_ratio\" : "
        << prior.aspect_ratio.value[0];
  }
  if (prior.skew.is_set) {
    out << ",\n";
    out << "\"skew\" : " << prior.skew.value[0];
  }
  if (prior.radial_distortion.is_set) {
    out << ",\n";
    out << "\"radial_distortion_coeffs\" : ["
        << prior.radial_distortion.value[0] << ", "
        << prior.radial_distortion.value[1] << ", "
        << prior.radial_distortion.value[2] << ", "
        << prior.radial_distortion.value[3] << "]";
  }
  if (prior.tangential_distortion.is_set) {
    out << ",\n";
    out << "\"tangential_distortion_coeffs\" : ["
        << prior.tangential_distortion.value[0] << ", "
        << prior.tangential_distortion.value[1] << "]";
  }
  if (prior.position.is_set) {
    out << ",\n";
    out << "\"position\" : ["
        << prior.position.value[0] << ", "
        << prior.position.value[1] << ", "
        << prior.position.value[2] << "]";
  }
  if (prior.orientation.is_set) {
    out << ",\n";
    out << "\"orientation\" : ["
        << prior.orientation.value[0] << ", "
        << prior.orientation.value[1] << ", "
        << prior.orientation.value[2] << "]";
  }
  if (prior.latitude.is_set) {
    out << ",\n";
    out << "\"latitude\" : " << prior.latitude.value[0];
  }
  if (prior.longitude.is_set) {
    out << ",\n";
    out << "\"longitude\" : " << prior.longitude.value[0];
  }
  if (prior.altitude.is_set) {
    out << ",\n";
    out << "\"altitude\" : " << prior.altitude.value[0];
  }
  out << "\n}}";
  return true;
}

}  // namespace

bool WriteCalibration(
    const std::string& output_calibration_file,
    const std::unordered_map<std::string, CameraIntrinsicsPrior>& priors) {
  // Return false if the unordered_map is empty.
  if (priors.empty()) {
    LOG(WARNING) << "The prior container is empty.";
    return false;
  }

  std::ofstream out(output_calibration_file);
  if (!out.is_open()) {
    LOG(ERROR) << "Could not write file: " << output_calibration_file;
    return false;
  }

  // Start the json file.
  out << "{\n\"priors\" : [\n";
  
  // Write each of the priors.
  std::unordered_map<std::string, CameraIntrinsicsPrior>::const_iterator it;
  it = priors.begin();
  // Write the first prior.
  WritePrior(it->first, it->second, &out);
  ++it;

  // Write the remaining priors.
  for ( ; it != priors.end(); ++it) {
    out << ",\n";
    WritePrior(it->first, it->second, &out);
  }

  // Close the json file.
  out << "\n]\n}\n";
  out.close();

  return true;
}

}  // namespace theia
