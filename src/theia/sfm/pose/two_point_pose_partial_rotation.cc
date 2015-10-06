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

#include "theia/sfm/pose/two_point_pose_partial_rotation.h"

#include <math.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include "theia/math/closed_form_polynomial_solver.h"

namespace theia {
using Eigen::Map;
using Eigen::Vector3d;
using Eigen::Quaterniond;

namespace {

// Adds a specific pose solution corresponding to the ray lengths.
inline void AddPoseSolution(const Vector3d& axis,
                            const Vector3d& model_point_1,
                            const Vector3d& model_point_2,
                            const Vector3d& image_ray_1,
                            const Vector3d& image_ray_2,
                            const double ray_length_1,
                            const double ray_length_2,
                            Quaterniond* rotation,
                            Vector3d* translation) {
  // Scales the image_rays by the calculated ray lengths, computing the points
  // in image space.
  const Vector3d point_in_image_space_1 = ray_length_1 * image_ray_1;
  const Vector3d point_in_image_space_2 = ray_length_2 * image_ray_2;
  const Vector3d image_points_diff =
      point_in_image_space_1 - point_in_image_space_2;

  const Vector3d model_points_diff = model_point_1 - model_point_2;

  // Computes two basis vectors that lie on the plane orthogonal to the axis.
  const Vector3d basis_vector_2 = axis.cross(model_points_diff).normalized();
  const Vector3d basis_vector_1 = basis_vector_2.cross(axis).normalized();

  DCHECK_LT(0.0, basis_vector_1.dot(basis_vector_1));

  // Finds the projection of the image_points_diff vector in this basis.
  const double dp_1 = basis_vector_1.dot(image_points_diff);
  const double dp_2 = basis_vector_2.dot(image_points_diff);

  // Finds the angle around the axis.
  const double angle = atan2(dp_2, dp_1);
  *rotation = Quaterniond(Eigen::AngleAxisd(angle, axis));

  // Calculates the translation, after taking into account the rotation.
  *translation = point_in_image_space_1 - *rotation * model_point_1;
}

int TwoPointPoseCore(const Vector3d& axis,
                     const Vector3d& model_point_1,
                     const Vector3d& model_point_2,
                     const Vector3d& image_ray_1,
                     const Vector3d& image_ray_2,
                     Quaterniond soln_rotations[2],
                     Vector3d soln_translations[2]) {
  // Let the points in the camera coordinate system be
  // y * image_ray_1, x * image_ray_2, where x, y are the lengths of the
  // rays. Since there is only rotation about the passed axis, the difference
  // between the values in the model and camera coordinate systems, projected
  // on the axis, should be the same. So:
  //
  // DotProd(axis, y * image_ray_1 - x * image_ray_2)  =
  //     DotProd(axis, model_points[0] - model_points[1])
  //
  // This allows x to be expressed in terms of y:
  //
  // x = m + n * y
  const double ray_1_axis_dp = image_ray_1.dot(axis);
  const double ray_2_axis_dp = image_ray_2.dot(axis);
  const double model_diff_axis_dp = (model_point_1 - model_point_2).dot(axis);

  const double m = model_diff_axis_dp / ray_1_axis_dp;
  const double n = ray_2_axis_dp / ray_1_axis_dp;

  // Next, the distance between the model points and the image points should
  // be the same:
  //
  // |y * image_ray_1 - x * image_ray_2| =
  //     |model_points[0] - model_points[1]|
  //
  // Using this and the substitution for x above we can create a quadratic
  // equation in y.

  // Computes the coefficients of the quadratic equation ay^2 + by + c = 0.0.
  const double ray_dp = image_ray_1.dot(image_ray_2);

  const long double a = n * (n - 2.0 * ray_dp) + 1.0;
  const long double b = 2.0 * m * (n - ray_dp);
  const long double c = m * m - (model_point_1 - model_point_2).squaredNorm();

  double roots[2] = { 0.0, 0.0 };
  const int number_of_roots = SolveQuadraticReals(a, b, c, roots);

  int num_solutions = 0;
  for (int i = 0; i < number_of_roots; ++i) {
    // Only accept positive ray distances.
    if (roots[i] > 0) {
      // Computes the other ray distance.
      const double ray_distance = m + n * roots[i];
      if (ray_distance > 0) {
        AddPoseSolution(axis, model_point_1, model_point_2, image_ray_1,
                        image_ray_2, ray_distance, roots[i],
                        &soln_rotations[num_solutions],
                        &soln_translations[num_solutions]);
        num_solutions++;
      }
    }
  }
  return num_solutions;
}

}  // namespace

int TwoPointPosePartialRotation(const Vector3d& axis,
                                const Vector3d& model_point_1,
                                const Vector3d& model_point_2,
                                const Vector3d& image_ray_1,
                                const Vector3d& image_ray_2,
                                Quaterniond soln_rotations[2],
                                Vector3d soln_translations[2]) {
  static const double kEpsilon = 1e-9;
  DCHECK_LT(fabs(image_ray_1.squaredNorm() - 1.0), kEpsilon);
  DCHECK_LT(fabs(image_ray_2.squaredNorm() - 1.0), kEpsilon);

  if (fabs(image_ray_1.dot(axis)) < kEpsilon) {
    if (fabs(image_ray_2.dot(axis)) > kEpsilon) {
      // If image_ray_1.y() == 0 then the function above doesn't work because
      // the calculate m and n will have a divide by 0.
      // However if image_ray_2.y() is not equal to zero then we can swap
      // and the points and call the above function with the swapped points.
      // TODO(cmsweeney): Maybe always swap to improve the conditioning?
      return TwoPointPoseCore(axis, model_point_2, model_point_1, image_ray_2,
                              image_ray_1, soln_rotations, soln_translations);
    }
    return 0;
  } else {
    return TwoPointPoseCore(axis, model_point_1, model_point_2, image_ray_1,
                            image_ray_2, soln_rotations, soln_translations);
  }
}

}  // namespace theia
