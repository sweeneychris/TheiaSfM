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

#include "theia/sfm/pose/dls_pnp.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <cmath>
#include <vector>

#include "theia/util/random.h"
#include "theia/sfm/pose/dls_impl.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector2d;
using Eigen::Vector3d;
using dls_impl::CreateMacaulayMatrix;
using dls_impl::ExtractJacobianCoefficients;
using dls_impl::LeftMultiplyMatrix;

// This implementation is ported from the Matlab version provided by the authors
// of "A Direct Least-Squares (DLS) Method for PnP". The general approach is to
// first rewrite the reprojection constraint (i.e., cost function) such that all
// unknowns appear linearly in terms of the rotation parameters (which are 3
// parameters in the Cayley-Gibss-Rodriguez formulation). Then we create a
// system of equations from the jacobian of the cost function, and solve these
// equations via a Macaulay matrix to obtain the roots (i.e., the 3 parameters
// of rotation). The translation can then be obtained through back-substitution.
void DlsPnp(const std::vector<Vector2d>& feature_position,
            const std::vector<Vector3d>& world_point,
            std::vector<Quaterniond>* solution_rotation,
            std::vector<Vector3d>* solution_translation) {
  CHECK_GE(feature_position.size(), 3);
  CHECK_EQ(feature_position.size(), world_point.size());

  const int num_correspondences = feature_position.size();

  // Holds the normalized feature positions cross multiplied with itself
  // i.e. n * n^t. This value is used multiple times so it is efficient to
  // pre-compute it.
  std::vector<Matrix3d> normalized_feature_cross;
  normalized_feature_cross.reserve(num_correspondences);
  for (int i = 0; i < num_correspondences; i++) {
    const Vector3d normalized_feature_pos =
        feature_position[i].homogeneous().normalized();
    normalized_feature_cross.push_back(normalized_feature_pos *
                                       normalized_feature_pos.transpose());
  }

  // The bottom-right symmetric block matrix of inverse(A^T * A). Matrix H from
  // Eq. 25 in the Appendix of the DLS paper.
  Matrix3d h_inverse = num_correspondences * Matrix3d::Identity();
  for (int i = 0; i < num_correspondences; i++) {
    h_inverse = h_inverse - normalized_feature_cross[i];
  }
  const Matrix3d h_matrix = h_inverse.inverse();

  // Compute V*W*b with the rotation parameters factored out. This is the
  // translation parameterized by the 9 entries of the rotation matrix.
  Matrix<double, 3, 9> translation_factor = Matrix<double, 3, 9>::Zero();
  for (int i = 0; i < num_correspondences; i++) {
    translation_factor = translation_factor +
                         (normalized_feature_cross[i] - Matrix3d::Identity()) *
                             LeftMultiplyMatrix(world_point[i]);
  }

  translation_factor = h_matrix * translation_factor;

  // Compute the cost function J' of Eq. 17 in DLS paper. This is a factorized
  // version where the rotation matrix parameters have been pulled out. The
  // entries to this equation are the coefficients to the cost function which is
  // a quartic in the rotation parameters.
  Matrix<double, 9, 9> ls_cost_coefficients = Matrix<double, 9, 9>::Zero();
  for (int i = 0; i < num_correspondences; i++) {
    ls_cost_coefficients =
        ls_cost_coefficients +
        (LeftMultiplyMatrix(world_point[i]) + translation_factor).transpose() *
            (Matrix3d::Identity() - normalized_feature_cross[i]) *
            (LeftMultiplyMatrix(world_point[i]) + translation_factor);
  }

  // Extract the coefficients of the jacobian (Eq. 18) from the
  // ls_cost_coefficients matrix. The jacobian represent 3 monomials in the
  // rotation parameters. Each entry of the jacobian will be 0 at the roots of
  // the polynomial, so we can arrange a system of polynomials from these
  // equations.
  double f1_coeff[20];
  double f2_coeff[20];
  double f3_coeff[20];
  ExtractJacobianCoefficients(ls_cost_coefficients, f1_coeff, f2_coeff,
                              f3_coeff);

  // We create one equation with random terms that is generally non-zero at the
  // roots of our system.
  const Eigen::Vector4d rand_vec = 100.0 * Eigen::Vector4d::Random();
  const double macaulay_term[4] = {rand_vec(0),
                                   rand_vec(1),
                                   rand_vec(2),
                                   rand_vec(3)};

  // Create Macaulay matrix that will be used to solve our polynonomial system.
  const MatrixXd& macaulay_matrix =
      CreateMacaulayMatrix(f1_coeff, f2_coeff, f3_coeff, macaulay_term);

  // Via the Schur complement trick, the top-left of the Macaulay matrix
  // contains a multiplication matrix whose eigenvectors correspond to solutions
  // to our system of equations.
  const MatrixXd solution_polynomial =
      macaulay_matrix.block<27, 27>(0, 0) -
      (macaulay_matrix.block<27, 93>(0, 27) *
       macaulay_matrix.block<93, 93>(27, 27).partialPivLu().solve(
           macaulay_matrix.block<93, 27>(27, 0)));

  // Extract eigenvectors of the solution polynomial to obtain the roots which
  // are contained in the entries of the eigenvectors.
  const Eigen::EigenSolver<MatrixXd> eigen_solver(solution_polynomial);

  // Many of the eigenvectors will contain complex solutions so we must filter
  // them to find the real solutions.
  const auto eigen_vectors = eigen_solver.eigenvectors();
  for (int i = 0; i < 27; i++) {
    // The first entry of the eigenvector should equal 1 according to our
    // polynomial, so we must divide each solution by the first entry.
    std::complex<double> s1 = eigen_vectors(9, i) / eigen_vectors(0, i);
    std::complex<double> s2 = eigen_vectors(3, i) / eigen_vectors(0, i);
    std::complex<double> s3 = eigen_vectors(1, i) / eigen_vectors(0, i);

    // If the rotation solutions are real, treat this as a valid candidate
    // rotation.
    const double kEpsilon = 1e-6;
    if (fabs(s1.imag()) < kEpsilon && fabs(s2.imag()) < kEpsilon &&
        fabs(s3.imag()) < kEpsilon) {
      // Compute the rotation (which is the transpose rotation of our solution)
      // and translation.
      Quaterniond soln_rotation(1.0, s1.real(), s2.real(), s3.real());
      soln_rotation = soln_rotation.inverse().normalized();

      const Matrix3d rot_mat = soln_rotation.inverse().toRotationMatrix();
      const Eigen::Map<const Matrix<double, 9, 1> > rot_vec(rot_mat.data());
      const Vector3d soln_translation = translation_factor * rot_vec;

      // TODO(cmsweeney): evaluate cost function and return it as an output
      // variable.

      // Check that all points are in front of the camera. Discard the solution
      // if this is not the case.
      bool all_points_in_front_of_camera = true;
      for (int j = 0; j < num_correspondences; j++) {
        const Vector3d transformed_point =
            soln_rotation * world_point[j] + soln_translation;
        if (transformed_point.z() < 0) {
          all_points_in_front_of_camera = false;
          break;
        }
      }

      if (all_points_in_front_of_camera) {
        solution_rotation->push_back(soln_rotation);
        solution_translation->push_back(soln_translation);
      }
    }
  }
}

}  // namespace theia
