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

#ifndef THEIA_SFM_CAMERA_INTRINSICS_PRIOR_H_
#define THEIA_SFM_CAMERA_INTRINSICS_PRIOR_H_

#include <cereal/access.hpp>

namespace theia {

// Weak calibration is not always available, so we need this helper struct to
// keep track of which data fields have been set.
struct Prior {
 public:
  bool is_set = false;
  double value = 0;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar) {  // NOLINT
    ar(is_set, value);
  }
};

// Prior information about a View. This is typically gathered from EXIF or
// sensor data that provides weak calibration.
struct CameraIntrinsicsPrior {
 public:
  // The image size *should* always be set, so we don't have to worry about
  // making it an Prior type.
  int image_width;
  int image_height;

  Prior focal_length;
  Prior principal_point[2];
  Prior aspect_ratio;
  Prior skew;
  Prior radial_distortion[2];

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar) {  // NOLINT
    ar(focal_length,
       principal_point[0],
       principal_point[1],
       aspect_ratio,
       skew,
       radial_distortion[0],
       radial_distortion[1]);
  }
};

}  // namespace theia

#endif  // THEIA_SFM_CAMERA_INTRINSICS_PRIOR_H_
