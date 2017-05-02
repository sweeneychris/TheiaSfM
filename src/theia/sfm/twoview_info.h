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

#ifndef THEIA_SFM_TWOVIEW_INFO_H_
#define THEIA_SFM_TWOVIEW_INFO_H_

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <stdint.h>
#include <Eigen/Core>

#include "theia/io/eigen_serializable.h"
#include "theia/sfm/types.h"

namespace theia {
class Camera;

// A struct to hold match and projection data between two views. It is assumed
// that the first view is at the origin with an identity rotation.
class TwoViewInfo {
 public:
  TwoViewInfo()
      : focal_length_1(0.0),
        focal_length_2(0.0),
        position_2(Eigen::Vector3d::Zero()),
        rotation_2(Eigen::Vector3d::Zero()),
        num_verified_matches(0),
        num_homography_inliers(0) {}

  double focal_length_1;
  double focal_length_2;

  Eigen::Vector3d position_2;
  Eigen::Vector3d rotation_2;

  // Number of features that were matched and geometrically verified betwen the
  // images.
  int num_verified_matches;

  // Number of inliers based on homography estimation. This is useful for
  // incremental SfM for choosing an initial view pair for the reconstruction.
  int num_homography_inliers;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {  // NOLINT
    ar(focal_length_1,
       focal_length_2,
       position_2,
       rotation_2,
       num_verified_matches,
       num_homography_inliers);
  }
};

// Inverts the two view info such that the focal lengths are swapped and the
// rotation and position are inverted.
void SwapCameras(TwoViewInfo* twoview_info);

// Constructs a TwoViewInfo object where camera1 is the "base" camera. In other
// words, the twoview info provides the relative pose of camera2 w.r.t. camera1.
void TwoViewInfoFromTwoCameras(const Camera& camera1,
                               const Camera& camera2,
                               TwoViewInfo* info);
}  // namespace theia

CEREAL_CLASS_VERSION(theia::TwoViewInfo, 0);

#endif  // THEIA_SFM_TWOVIEW_INFO_H_
