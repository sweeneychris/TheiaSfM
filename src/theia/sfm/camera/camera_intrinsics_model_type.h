// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_CAMERA_CAMERA_INTRINSICS_MODEL_TYPE_H_
#define THEIA_SFM_CAMERA_CAMERA_INTRINSICS_MODEL_TYPE_H_

#include <glog/logging.h>
#include <string>

namespace theia {

// Each camera model implemented through this interface should have a type
// listed here. The Create method below should create an instance to the
// respective camera model based on the type provided.
enum class CameraIntrinsicsModelType {
  INVALID = -1,
  PINHOLE = 0,
  PINHOLE_RADIAL_TANGENTIAL = 1,
  FISHEYE = 2,
  FOV = 3,
  DIVISION_UNDISTORTION = 4,
};

// Converts an input string to the corresponding camera intrinsics model type.
inline CameraIntrinsicsModelType StringToCameraIntrinsicsModelType(
    const std::string& camera_model_type_string) {
  if (camera_model_type_string == "PINHOLE") {
    return CameraIntrinsicsModelType::PINHOLE;
  } else if (camera_model_type_string == "PINHOLE_RADIAL_TANGENTIAL") {
    return CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL;
  } else if (camera_model_type_string == "FISHEYE") {
    return CameraIntrinsicsModelType::FISHEYE;
  } else if (camera_model_type_string == "FOV") {
    return CameraIntrinsicsModelType::FOV;
  } else if (camera_model_type_string == "DIVISION_UNDISTORTION") {
    return CameraIntrinsicsModelType::DIVISION_UNDISTORTION;
  } else {
    LOG(FATAL) << "Invalid camera model type supplied: "
               << camera_model_type_string;
  }
}

inline std::string CameraIntrinsicsModelTypeToString(
    const CameraIntrinsicsModelType& camera_model_type) {
  switch (camera_model_type) {
    case CameraIntrinsicsModelType::PINHOLE:
      return "PINHOLE";
    case CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL:
      return "PINHOLE_RADIAL_TANGENTIAL";
    case CameraIntrinsicsModelType::FISHEYE:
      return "FISHEYE";
    case CameraIntrinsicsModelType::FOV:
      return "FOV";
    case CameraIntrinsicsModelType::DIVISION_UNDISTORTION:
      return "DIVISION_UNDISTORTION";
    default:
      LOG(FATAL) << "Invalid Camera model chosen.";
      break;
  }

  LOG(FATAL) << "This should not be reached!!";
  return "";
}

// Converts an input string to the corresponding camera intrinsics model type.
inline bool IsCameraIntrinsicsModelTypeValid(
    const std::string& camera_model_type_string) {
  if (camera_model_type_string == "PINHOLE") {
    return true;
  } else if (camera_model_type_string == "PINHOLE_RADIAL_TANGENTIAL") {
    return true;
  } else if (camera_model_type_string == "FISHEYE") {
    return true;
  } else if (camera_model_type_string == "FOV") {
    return true;
  } else if (camera_model_type_string == "DIVISION_UNDISTORTION") {
    return true;
  }
  return false;
}

}  // namespace theia

#endif  // THEIA_SFM_CAMERA_CAMERA_INTRINSICS_MODEL_TYPE_H_
