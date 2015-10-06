// Copyright (C) 2015 The Regents of the University of California (Regents)
// and Google, Inc. All rights reserved.
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
//     * Neither the name of The Regents or University of California, Google,
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
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

#ifndef THEIA_SFM_POSE_TWO_POINT_POSE_PARTIAL_ROTATION_H_
#define THEIA_SFM_POSE_TWO_POINT_POSE_PARTIAL_ROTATION_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace theia {
// Solves for the limited pose of a camera from two 3D points to image ray
// correspondences. The pose is limited in that while it solves for the three
// translation components, it only solves for a single rotation around a passed
// axis.
//
// This is intended for use with camera phones that have accelerometers, so that
// the 'up' vector is known, meaning the other two rotations are known. The
// effect of the other rotations should be removed before using this function.
//
// This implementation is intended to form the core of a RANSAC routine, and as
// such has an optimized interface for this use case.
//
// Computes the limited pose between the 3D model points and the (unit-norm)
// image rays. Places the rotation and translation solutions in soln_rotations
// and soln_translations.
// There are at most 2 solutions, and the number of solutions is returned.
//
// The rotations and translation are defined such that model points are
// transformed according to:
//
//     image_point = Q * model_point + t
//
// This function computes the rotation and translation such that the model
// points, after transformation, lie along the corresponding image_rays. The
// axis referred to is the axis of rotation between the camera coordinate system
// and world (3D point) coordinate system. For most users, this axis will be
// (0, 1, 0) i.e., the up direction. This requires that the input image rays
// have been rotated such that the up direction of the camera coordinate system
// is indeed equal to (0, 1, 0).
//
// When using this algorithm please cite the paper "Efficient Computation of
// Absolute Pose for Gravity-Aware Augmented Reality" by Sweeney et al (ISMAR
// 2015).
int TwoPointPosePartialRotation(const Eigen::Vector3d& axis,
                                const Eigen::Vector3d& model_point_1,
                                const Eigen::Vector3d& model_point_2,
                                const Eigen::Vector3d& image_ray_1,
                                const Eigen::Vector3d& image_ray_2,
                                Eigen::Quaterniond soln_rotations[2],
                                Eigen::Vector3d soln_translations[2]);
}  // namespace theia

#endif  // THEIA_SFM_POSE_TWO_POINT_POSE_PARTIAL_ROTATION_H_
