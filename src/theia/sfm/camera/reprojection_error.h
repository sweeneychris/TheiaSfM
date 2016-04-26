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

#ifndef THEIA_SFM_CAMERA_REPROJECTION_ERROR_H_
#define THEIA_SFM_CAMERA_REPROJECTION_ERROR_H_

#include <ceres/ceres.h>
#include "theia/sfm/feature.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/project_point_to_image.h"

namespace theia {

struct ReprojectionError {
 public:
  explicit ReprojectionError(const Feature& feature) : feature_(feature) {}

  template<typename T> bool operator()(const T* camera_extrinsics,
                                       const T* camera_intrinsics,
                                       const T* point_parameters,
                                       T* reprojection_error) const {
    // Do not evaluate invalid camera configurations.
    if (camera_intrinsics[Camera::FOCAL_LENGTH] < T(0.0) ||
        camera_intrinsics[Camera::ASPECT_RATIO] < T(0.0)) {
      return false;
    }

    T reprojection[2];
    ProjectPointToImage(camera_extrinsics,
                        camera_intrinsics,
                        point_parameters,
                        reprojection);
    reprojection_error[0] = reprojection[0] - T(feature_.x());
    reprojection_error[1] = reprojection[1] - T(feature_.y());
    return true;
  }

  static ceres::CostFunction* Create(const Feature& feature) {
    static const int kPointSize = 4;
    return new ceres::AutoDiffCostFunction<ReprojectionError,
                                           2,
                                           Camera::kExtrinsicsSize,
                                           Camera::kIntrinsicsSize,
                                           kPointSize>(
        new ReprojectionError(feature));
  }

 private:
  const Feature feature_;
};

}  // namespace theia

#endif  // THEIA_SFM_CAMERA_REPROJECTION_ERROR_H_
