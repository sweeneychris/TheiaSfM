// Copyright (C) 2018 The Regents of the University of California (Regents).
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

#include "theia/sfm/pose/upnp.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include <glog/logging.h>

#include "theia/alignment/alignment.h"

namespace theia {

namespace {
typedef Eigen::Matrix<double, 3, 10> Matrix3x10d;

// Computes the H Matrix (see Eq. (6)) and the outer products of the ray
// directions, since these are used to compute matrix V (Eq. (5)).
inline Eigen::Matrix3d ComputeHMatrixAndRayDirectionsOuterProducts(
    const std::vector<Eigen::Vector3d>& ray_directions,
    std::vector<Eigen::Matrix3d>* outer_products) {
  CHECK_NOTNULL(outer_products)->reserve(ray_directions.size());
  Eigen::Matrix3d h_inverse;
  h_inverse.setZero();
  for (const Eigen::Vector3d& ray : ray_directions) {
    outer_products->emplace_back(ray * ray.transpose());
    h_inverse -= outer_products->back();
  }
  h_inverse += ray_directions.size() * Eigen::Matrix3d::Identity();
  return h_inverse.inverse();
}

void LeftMultiply(const Eigen::Vector3d& point,
                  Matrix3x10d* left_multiply_matrix) {
  Matrix3x10d& phi_mat = *CHECK_NOTNULL(left_multiply_matrix);
  // Row 0.
  phi_mat(0, 0) = point.x();
  phi_mat(0, 1) = point.x();
  phi_mat(0, 2) = -point.x();
  phi_mat(0, 3) = -point.x();
  phi_mat(0, 4) = 0.0;
  phi_mat(0, 5) = 2 * point.z();
  phi_mat(0, 6) = -2 * point.y();
  phi_mat(0, 7) = 2 * point.y();
  phi_mat(0, 8) = 2 * point.z();
  phi_mat(0, 9) = 0.0;

  // Row 1.
  phi_mat(1, 0) = point.y();
  phi_mat(1, 1) = -point.y();
  phi_mat(1, 2) = point.y();
  phi_mat(1, 3) = -point.y();
  phi_mat(1, 4) = -2.0 * point.z();
  phi_mat(1, 5) = 0.0;
  phi_mat(1, 6) = 2 * point.x();
  phi_mat(1, 7) = 2 * point.x();
  phi_mat(1, 8) = 0.0;
  phi_mat(1, 9) = 2 * point.z();

  // Row 3.
  phi_mat(2, 0) = point.y();
  phi_mat(2, 1) = -point.y();
  phi_mat(2, 2) = point.y();
  phi_mat(2, 3) = -point.y();
  phi_mat(2, 4) = -2.0 * point.z();
  phi_mat(2, 5) = 0.0;
  phi_mat(2, 6) = 2 * point.x();
  phi_mat(2, 7) = 2 * point.x();
  phi_mat(2, 8) = 0.0;
  phi_mat(2, 9) = 2 * point.z();
}

void ComputeHelperMatrices(
    const std::vector<Eigen::Vector3d>& world_points,
    const std::vector<Eigen::Vector3d>& ray_origins,
    const std::vector<Eigen::Matrix3d>& outer_products,
    const Eigen::Matrix3d& h_matrix,
    Matrix3x10d* g_matrix,
    Eigen::Vector3d* j_matrix) {
  CHECK_EQ(ray_origins.size(), outer_products.size());
  CHECK_NOTNULL(g_matrix)->setZero();
  CHECK_NOTNULL(j_matrix)->setZero();
  const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
  Matrix3x10d left_multiply_mat;
  for (int i = 0; i < ray_origins.size(); ++i) {
    const Eigen::Matrix3d& outer_product = outer_products[i];
    // Computation following Eq. (5).
    const Eigen::Matrix3d v_matrix = h_matrix * (outer_product - identity);
    // Compute the left multiplication matrix or Phi matrix in the paper.
    LeftMultiply(world_points[i], &left_multiply_mat);
    *j_matrix += v_matrix * ray_origins[i];
    *g_matrix += v_matrix * left_multiply_mat;
  }
}

}  // namespace

// TODO(vfragoso): Document me!
void Upnp(const std::vector<Eigen::Vector3d>& ray_origins,
          const std::vector<Eigen::Vector3d>& ray_directions,
          const std::vector<Eigen::Vector3d>& world_points,
          std::vector<Eigen::Quaterniond>* solution_rotations,
          std::vector<Eigen::Vector3d>* solution_translations) {
  CHECK_NOTNULL(solution_rotations)->clear();
  CHECK_NOTNULL(solution_translations)->clear();

  // 1. Compute the H matrix and the outer products of the ray directions.
  std::vector<Eigen::Matrix3d> outer_products;
  const Eigen::Matrix3d h_matrix =
      ComputeHMatrixAndRayDirectionsOuterProducts(
          ray_directions, &outer_products);

  // 2. Compute matrices J and G from page 132 or 6-th page in the paper.
  Matrix3x10d g_matrix;
  Eigen::Vector3d j_matrix;
  ComputeHelperMatrices(world_points,
                        ray_origins,
                        outer_products,
                        h_matrix,
                        &g_matrix,
                        &j_matrix);
}

}  // namespace theia

