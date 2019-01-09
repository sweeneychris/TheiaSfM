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
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include <complex>
#include <vector>

#include "theia/alignment/alignment.h"
#include "theia/sfm/pose/build_upnp_action_matrix.h"

namespace theia {

namespace {
typedef Eigen::Matrix<double, 3, 10> Matrix3x10d;
typedef Eigen::Matrix<double, 8, 8> Matrix8d;
typedef Eigen::Matrix<std::complex<double>, 8, 8> Matrix8cd;
typedef Eigen::Matrix<double, 10, 10> Matrix10d;
typedef Eigen::Matrix<double, 16, 16> Matrix16d;
typedef Eigen::Matrix<std::complex<double>, 16, 16> Matrix16cd;
typedef Eigen::Matrix<double, 10, 1> Vector10d;

const int kNumMaxRotations = 16;
const int kNumMaxRotationsExploitingSymmetry = 8;
const int kNumMinCorrespondences = 4;

// TODO(vfragoso): Document me!
struct InputDatum {
  InputDatum(const std::vector<Eigen::Vector3d>& _ray_origins,
             const std::vector<Eigen::Vector3d>& _ray_directions,
             const std::vector<Eigen::Vector3d>& _world_points) :
      ray_origins(_ray_origins),
      ray_directions(_ray_directions),
      world_points(_world_points) {}
  ~InputDatum() = default;

  const std::vector<Eigen::Vector3d>& ray_origins;
  const std::vector<Eigen::Vector3d>& ray_directions;
  const std::vector<Eigen::Vector3d>& world_points;
};

// Computes the H Matrix (see Eq. (6)) and the outer products of the ray
// directions, since these are used to compute matrix V (Eq. (5)).
inline Eigen::Matrix3d ComputeHMatrixAndRayDirectionsOuterProducts(
    const InputDatum& input_datum,
    std::vector<Eigen::Matrix3d>* outer_products) {
  const std::vector<Eigen::Vector3d>& ray_directions =
      input_datum.ray_directions;
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

inline Matrix3x10d LeftMultiply(const Eigen::Vector3d& point) {
  Matrix3x10d phi_mat;
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
  phi_mat(2, 0) = point.z();
  phi_mat(2, 1) = -point.z();
  phi_mat(2, 2) = -point.z();
  phi_mat(2, 3) = point.z();
  phi_mat(2, 4) = 2.0 * point.y();
  phi_mat(2, 5) = -2.0 * point.x();
  phi_mat(2, 6) = 0.0;
  phi_mat(2, 7) = 0.0;
  phi_mat(2, 8) = 2.0 * point.x();
  phi_mat(2, 9) = 2.0 * point.y();
  return phi_mat;
}

inline void ComputeHelperMatrices(
    const InputDatum& input_datum,
    const std::vector<Eigen::Matrix3d>& outer_products,
    const Eigen::Matrix3d& h_matrix,
    Matrix3x10d* g_matrix,
    Eigen::Vector3d* j_matrix) {
  const std::vector<Eigen::Vector3d>& world_points = input_datum.world_points;
  const std::vector<Eigen::Vector3d>& ray_origins = input_datum.ray_origins;
  CHECK_EQ(ray_origins.size(), outer_products.size());
  CHECK_NOTNULL(g_matrix)->setZero();
  CHECK_NOTNULL(j_matrix)->setZero();
  const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
  for (int i = 0; i < ray_origins.size(); ++i) {
    const Eigen::Matrix3d& outer_product = outer_products[i];
    // Computation following Eq. (5).
    const Eigen::Matrix3d v_matrix = h_matrix * (outer_product - identity);
    // Compute the left multiplication matrix or Phi matrix in the paper.
    const Matrix3x10d left_multiply_mat = LeftMultiply(world_points[i]);
    *j_matrix += v_matrix * ray_origins[i];
    *g_matrix += v_matrix * left_multiply_mat;
  }
}

// Computes the block matrices that compose the M matrix in Eq. 17. These
// blocks are:
// a_matrix = \sum A_i^T * A_i,
// b_vector = \sum A_i^T * b_i ,
// gamma = \sum b_i^T * b_i.
UpnpCostParameters ComputeCostParameters(
    const InputDatum& input_datum,
    const std::vector<Eigen::Matrix3d>& outer_products,
    const Matrix3x10d& g_matrix,
    const Eigen::Vector3d& j_matrix) {
  const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
  UpnpCostParameters cost_params;
  const std::vector<Eigen::Vector3d>& world_points = input_datum.world_points;
  const std::vector<Eigen::Vector3d>& ray_origins = input_datum.ray_origins;
  Matrix10d& a_matrix = cost_params.a_matrix;
  Vector10d& b_vector = cost_params.b_vector;
  double& gamma = cost_params.gamma;
  // Gamma is the sum of the dot products of b_matrices.
  for (int i = 0; i < world_points.size(); ++i) {
    // Compute the left multiplication matrix or Phi matrix in the paper.
    const Matrix3x10d left_multiply_mat = LeftMultiply(world_points[i]);
    const Eigen::Matrix3d outer_prod_minus_identity =
        outer_products[i] - identity;
    // Compute the i-th a_matrix.
    const Matrix3x10d temp_a_mat =
        outer_prod_minus_identity * (left_multiply_mat + g_matrix);
    a_matrix += temp_a_mat.transpose() * temp_a_mat;
    // Compute the i-th b_vector.
    const Eigen::Vector3d temp_b_mat =
        -outer_prod_minus_identity * (ray_origins[i] + j_matrix);
    b_vector += temp_a_mat.transpose() * temp_b_mat;
    // Compute the i-th gamma.
    gamma += temp_b_mat.squaredNorm();
  }

  return cost_params;
}

std::vector<Eigen::Quaterniond> SolveUpnpFromNonMinimalSample(
    const InputDatum& input_datum,
    const UpnpCostParameters& cost_params) {
  std::vector<Eigen::Quaterniond> rotations(kNumMaxRotationsExploitingSymmetry);
  // Build action matrix.
  const Matrix8d action_matrix = BuildActionMatrixUsingSymmetry(
      cost_params.a_matrix,
      cost_params.b_vector,
      cost_params.gamma);
  const Eigen::EigenSolver<Matrix8d> eigen_solver(action_matrix, true);
  const Matrix8cd eigen_vectors = eigen_solver.eigenvectors();
  for (int i = 0; i < rotations.size(); ++i) {
    Eigen::Quaterniond quaternion(eigen_vectors(4, i).real(),
                                  eigen_vectors(5, i).real(),
                                  eigen_vectors(6, i).real(),
                                  eigen_vectors(7, i).real());
    quaternion.normalize();
    rotations[i] = quaternion;
  }

  return rotations;
}

std::vector<Eigen::Quaterniond> SolveUpnpFromMinimalSample(
    const InputDatum& input_datum,
    const UpnpCostParameters& cost_params) {
  std::vector<Eigen::Quaterniond> rotations(kNumMaxRotations);
  // Build action matrix.
  const Matrix16d action_matrix = BuildActionMatrix(cost_params.a_matrix,
                                                    cost_params.b_vector,
                                                    cost_params.gamma);
  const Eigen::EigenSolver<Matrix16d> eigen_solver(action_matrix, true);
  const Matrix16cd eigen_vectors = eigen_solver.eigenvectors();
  for (int i = 0; i < rotations.size(); ++i) {
    // According to the original implementation, the complex solutions
    // can be good, in particular when the number of correspondences is really
    // low. The solutions simply ignore the imaginary part.
    Eigen::Vector4d quaternion(eigen_vectors(11, i).real(),
                               eigen_vectors(12, i).real(),
                               eigen_vectors(13, i).real(),
                               eigen_vectors(14, i).real());
    quaternion.normalize();

    if (quaternion[0] < 0.0) {
      quaternion *= -1.0;
    }

    rotations[i] = Eigen::Quaterniond(quaternion[0],
                                      quaternion[1],
                                      quaternion[2],
                                      quaternion[3]);
  }

  return rotations;
}

inline std::vector<Eigen::Quaterniond> ComputeRotations(
    const InputDatum& input_datum,
    const UpnpCostParameters& cost_params) {
  // Build the action matrix.
  if (input_datum.world_points.size() > kNumMinCorrespondences) {
    return SolveUpnpFromNonMinimalSample(input_datum, cost_params);
  }
  return SolveUpnpFromMinimalSample(input_datum, cost_params);
}

// Constructs the vector s as indicated in Eq. 12.
inline Vector10d ComputeRotationVector(const Eigen::Quaterniond& rotation) {
  Vector10d rotation_vector;
  // Set the values of the rotation vector.
  rotation_vector[0] = rotation.w() * rotation.w();
  rotation_vector[1] = rotation.x() * rotation.x();
  rotation_vector[2] = rotation.y() * rotation.y();
  rotation_vector[3] = rotation.z() * rotation.z();
  rotation_vector[4] = rotation.w() * rotation.x();
  rotation_vector[5] = rotation.w() * rotation.y();
  rotation_vector[6] = rotation.w() * rotation.z();
  rotation_vector[7] = rotation.x() * rotation.y();
  rotation_vector[8] = rotation.x() * rotation.z();
  rotation_vector[9] = rotation.y() * rotation.z();
  return rotation_vector;
}

}  // namespace

// Evaluates the cost introduced in Eq. 17.
double EvaluateUpnpCost(const UpnpCostParameters& parameters,
                        const Eigen::Quaterniond& rotation) {
  // Compute the quaternion vector.
  const Vector10d rotation_vector = ComputeRotationVector(rotation);
  const double cost =
      (rotation_vector.transpose() * parameters.a_matrix * rotation_vector +
       2.0 * parameters.b_vector.transpose() * rotation_vector)(0, 0) +
      parameters.gamma;
  return cost;
}

UpnpCostParameters ComputeUpnpCostParameters(
    const std::vector<Eigen::Vector3d>& ray_origins,
    const std::vector<Eigen::Vector3d>& ray_directions,
    const std::vector<Eigen::Vector3d>& world_points) {
  const InputDatum input_datum(ray_origins, ray_directions, world_points);
  // 1. Compute the H matrix and the outer products of the ray directions.
  std::vector<Eigen::Matrix3d> outer_products;
  const Eigen::Matrix3d h_matrix =
      ComputeHMatrixAndRayDirectionsOuterProducts(input_datum, &outer_products);

  // 2. Compute matrices J and G from page 132 or 6-th page in the paper.
  Matrix3x10d g_matrix;
  Eigen::Vector3d j_matrix;
  ComputeHelperMatrices(input_datum,
                        outer_products,
                        h_matrix,
                        &g_matrix,
                        &j_matrix);

  // 3. Compute matrix the block-matrix of matrix M from Eq. 17.
  return ComputeCostParameters(input_datum,
                               outer_products,
                               g_matrix,
                               j_matrix);
}

// TODO(vfragoso): Document me!
UpnpCostParameters Upnp(const std::vector<Eigen::Vector3d>& ray_origins,
                        const std::vector<Eigen::Vector3d>& ray_directions,
                        const std::vector<Eigen::Vector3d>& world_points,
                        std::vector<Eigen::Quaterniond>* solution_rotations,
                        std::vector<Eigen::Vector3d>* solution_translations) {
  CHECK_NOTNULL(solution_rotations)->clear();
  CHECK_NOTNULL(solution_translations)->clear();

  // Compute Upnp cost parameters.
  const InputDatum input_datum(ray_origins, ray_directions, world_points);
  const UpnpCostParameters cost_params =
      ComputeUpnpCostParameters(ray_origins, ray_directions, world_points);

  // Compute rotations.
  const std::vector<Eigen::Quaterniond> rotations =
      ComputeRotations(input_datum, cost_params);

  // TODO(vfragoso): Compute translation.

  *solution_rotations = std::move(rotations);

  return cost_params;
}

}  // namespace theia

