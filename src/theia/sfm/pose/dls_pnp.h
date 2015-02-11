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

#ifndef THEIA_SFM_POSE_DLS_PNP_H_
#define THEIA_SFM_POSE_DLS_PNP_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "theia/alignment/alignment.h"

namespace theia {
// Computes the camera pose using the Perspective N-point method from "A Direct
// Least-Squares (DLS) Method for PnP" by Joel Hesch and Stergios
// Roumeliotis. This method is extremely scalable and highly accurate for the
// PnP problem. A minimum of 3 points are required, but there is no maximum
// number of points allowed as this is a least-squared approach. Theoretically,
// up to 27 solutions may be returned, but in practice only 4 real solutions
// arise and in almost all cases where n >= 6 there is only one solution which
// places the observed points in front of the camera.
//
// Params:
//   feature_position: Feature positions corresponding to model points. Must
//     contain at least 3 points.
//   points_3d: 3D location of features. Must correspond to the image_ray
//     of the same index. Must contain the same number of points as image_ray,
//     and at least 3.
//   solution_rotation: the rotation quaternion of the candidate solutions
//   solution_translation: the translation of the candidate solutions
void DlsPnp(const std::vector<Eigen::Vector2d>& feature_positions,
            const std::vector<Eigen::Vector3d>& world_point,
            std::vector<Eigen::Quaterniond>* solution_rotation,
            std::vector<Eigen::Vector3d>* solution_translation);
}  // namespace theia

#endif  // THEIA_SFM_POSE_DLS_PNP_H_
