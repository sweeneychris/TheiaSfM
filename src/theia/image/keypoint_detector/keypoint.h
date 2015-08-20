// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#ifndef THEIA_IMAGE_KEYPOINT_DETECTOR_KEYPOINT_H_
#define THEIA_IMAGE_KEYPOINT_DETECTOR_KEYPOINT_H_

#include <vector>

namespace theia {
// A generic keypoint class that mimics a protocol buffer The only variable
// thats must be set are x, y, and the type. All other variable are optional and
// will be set to THEIA_INVALID_KEYPOINT_VAR by default.
#define THEIA_INVALID_KEYPOINT_VAR -9999
class Keypoint {
 public:
  enum KeypointType {
    INVALID = -1,
    OTHER = 0,
    SIFT,
  };

  Keypoint(double x, double y, KeypointType type)
      : x_(x), y_(y), keypoint_type_(type),
        strength_(THEIA_INVALID_KEYPOINT_VAR),
        scale_(THEIA_INVALID_KEYPOINT_VAR),
        orientation_(THEIA_INVALID_KEYPOINT_VAR) {}

  Keypoint() : Keypoint(THEIA_INVALID_KEYPOINT_VAR,
                        THEIA_INVALID_KEYPOINT_VAR,
                        Keypoint::INVALID) {}

  ~Keypoint() {}

  // Keypoint type.
  inline KeypointType keypoint_type() const { return keypoint_type_; }
  inline void set_keypoint_type(KeypointType type) { keypoint_type_ = type; }

  // Variable x.
  inline double x() const { return x_; }
  inline void set_x(double x) { x_ = x; }

  // Variable y.
  inline double y() const { return y_; }
  inline void set_y(double y) { y_ = y; }

  // Optional variable strength.
  inline bool has_strength() const {
    return strength_ != THEIA_INVALID_KEYPOINT_VAR; }
  inline double strength() const { return strength_; }
  inline void set_strength(double strength) { strength_ = strength; }

  // Optional variable scale.
  inline bool has_scale() const { return scale_ != THEIA_INVALID_KEYPOINT_VAR; }
  inline double scale() const { return scale_; }
  inline void set_scale(double scale) { scale_ = scale; }

  // Optional variable orientation.
  inline bool has_orientation() const {
    return orientation_ != THEIA_INVALID_KEYPOINT_VAR; }
  inline double orientation() const { return orientation_; }
  inline void set_orientation(double orientation) {
    orientation_ = orientation; }

 private:
  double x_;
  double y_;
  KeypointType keypoint_type_;
  double strength_;
  double scale_;
  double orientation_;
};

}  // namespace theia

#endif  // THEIA_IMAGE_KEYPOINT_DETECTOR_KEYPOINT_H_
