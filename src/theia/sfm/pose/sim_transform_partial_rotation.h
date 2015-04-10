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

#ifndef THEIA_SFM_POSE_SIM_TRANSFORM_PARTIAL_ROTATION_H_
#define THEIA_SFM_POSE_SIM_TRANSFORM_PARTIAL_ROTATION_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>

namespace theia {

// Solves for the similarity transformation that will transform rays in image
// two such that the intersect with rays in image one such that:
//   s * R * X' + t = X
// where s, R, t are the scale, rotation, and translation returned, X' is a
// point in coordinate system 2 and X is the point transformed back to
// coordinate system 1.
//
// Please cite the paper "Computing Similarity Transformations from Only Image
// Correspondences" by C. Sweeney et al (CVPR 2015) when using this algorithm.
void SimTransformPartialRotation(
    const Eigen::Vector3d& rotation_axis,
    const Eigen::Vector3d image_one_ray_directions[5],
    const Eigen::Vector3d image_one_ray_origins[5],
    const Eigen::Vector3d image_two_ray_directions[5],
    const Eigen::Vector3d image_two_ray_origins[5],
    std::vector<Eigen::Quaterniond>* soln_rotations,
    std::vector<Eigen::Vector3d>* soln_translations,
    std::vector<double>* soln_scales);

}  // namespace theia

#endif  // THEIA_SFM_POSE_SIM_TRANSFORM_PARTIAL_ROTATION_H_
