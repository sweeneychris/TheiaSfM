// Copyright (C) 2017 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#include "theia/math/constrained_l1_solver.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <glog/logging.h>

#include <algorithm>
#include <string>

#include "theia/math/matrix/sparse_cholesky_llt.h"
#include "theia/util/stringprintf.h"

namespace theia {

ConstrainedL1Solver::ConstrainedL1Solver(
    const Options& options,
    const Eigen::SparseMatrix<double>& A,
    const Eigen::VectorXd& b,
    const Eigen::SparseMatrix<double>& geq_mat,
    const Eigen::VectorXd& geq_vec)
    : options_(options),
      num_l1_residuals_(b.size()),
      num_inequality_constraints_(geq_vec.size()) {
  CHECK_EQ(A.cols(), geq_mat.cols());
  CHECK_EQ(A.rows(), b.rows());
  CHECK_EQ(geq_mat.rows(), geq_vec.rows());

  // Allocate matrix A.
  A_.resize(A.rows() + geq_mat.rows(), A.cols());

  // Iterate over the input mat and geq_mat and store the entries in A.
  std::vector<Eigen::Triplet<double> > triplets;
  triplets.reserve(A.nonZeros() + geq_mat.nonZeros());
  for (int i = 0; i < A.outerSize(); i++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
      triplets.emplace_back(it.row(), it.col(), it.value());
    }
  }
  // Store the Geq mat below matrix A.
  for (int i = 0; i < geq_mat.outerSize(); i++) {
    // Iterate over inside
    for (Eigen::SparseMatrix<double>::InnerIterator it(geq_mat, i); it; ++it) {
      triplets.emplace_back(A.rows() + it.row(), it.col(), it.value());
    }
  }
  A_.setFromTriplets(triplets.begin(), triplets.end());

  Eigen::SparseMatrix<double> spd_mat(A.cols(), A.cols());
  spd_mat.selfadjointView<Eigen::Upper>().rankUpdate(A_.transpose());

  linear_solver_.Compute(spd_mat);
  CHECK_EQ(linear_solver_.Info(), Eigen::Success);

  // Set the modified b vector.
  b_.resize(b.size() + geq_vec.size());
  b_.head(b.size()) = b;
  b_.tail(geq_vec.size()) = geq_vec;
}

// We create a modified L1 solver such that ||Bx - b|| is minimized under L1
// norm subject to the constraint geq_mat * x > geq_vec. We conveniently
// create this constraint in ADMM terms as:
//
//    minimize f(x) + g(z_1) + h(z_2)
//    s.t. Bx - b - z_1 = 0
//         Cx - c - z_2 = 0
//
// Where f(x) = 0, g(z_1) = |z_1| and h(z_2) is an indicate function for our
// inequality constraint. This can be transformed into the standard ADMM
// formulation as:
//
//    minimize f(x) + g(z)
//    s.t. A * x - d - z = 0
//
// where A = [B;C] and d=[b;c] (where ; is the "stack" operation like matlab)
// This can now be solved in the same form as the L1 minimization, with a
// slightly different z update.
void ConstrainedL1Solver::Solve(Eigen::VectorXd* solution) {
  CHECK_NOTNULL(solution)->resize(A_.cols());
  Eigen::VectorXd& x = *solution;
  Eigen::VectorXd z(A_.rows()), u(A_.rows());
  z.setZero();
  u.setZero();

  Eigen::VectorXd a_times_x(A_.rows()), z_old(z.size()), ax_hat(A_.rows());
  // Precompute some convergence terms.
  const double rhs_norm = b_.norm();
  const double primal_abs_tolerance_eps =
      std::sqrt(A_.rows()) * options_.absolute_tolerance;
  const double dual_abs_tolerance_eps =
      std::sqrt(A_.cols()) * options_.absolute_tolerance;
  VLOG(2) << "Iteration   R norm          S norm          Primal eps      "
             "Dual eps";
  const std::string row_format =
      "  % 4d     % 4.4e     % 4.4e     % 4.4e     % 4.4e";

  // qp_options.max_num_iterations = 100;
  for (int i = 0; i < options_.max_num_iterations; i++) {
    x.noalias() = linear_solver_.Solve(A_.transpose() * (b_ + z - u));

    if (linear_solver_.Info() != Eigen::Success) {
      LOG(ERROR) << "L1 Minimization failed. Could not solve the sparse "
                    "linear system with Cholesky Decomposition";
      return;
    }

    a_times_x.noalias() = A_ * x;
    ax_hat.noalias() = options_.alpha * a_times_x;
    ax_hat.noalias() += (1.0 - options_.alpha) * (z + b_);

    // Update z and set z_old.
    std::swap(z, z_old);
    z.noalias() = ModifiedShrinkage(ax_hat - b_ + u, 1.0 / options_.rho);

    // Update u.
    u.noalias() += ax_hat - z - b_;

    // Compute the convergence terms.
    const double r_norm = (a_times_x - z - b_).norm();
    const double s_norm = (-options_.rho * A_.transpose() * (z - z_old)).norm();
    const double max_norm = std::max({a_times_x.norm(), z.norm(), rhs_norm});
    const double primal_eps =
        primal_abs_tolerance_eps + options_.relative_tolerance * max_norm;
    const double dual_eps = dual_abs_tolerance_eps +
                            options_.relative_tolerance *
                                (options_.rho * A_.transpose() * u).norm();

    // Log the result to the screen.
    VLOG(2) << theia::StringPrintf(
        row_format.c_str(), i, r_norm, s_norm, primal_eps, dual_eps);
    // Determine if the minimizer has converged.
    if (r_norm < primal_eps && s_norm < dual_eps) {
      break;
    }
  }
}

Eigen::VectorXd ConstrainedL1Solver::ModifiedShrinkage(
    const Eigen::VectorXd& vec, const double kappa) {
  Eigen::VectorXd output(num_l1_residuals_ + num_inequality_constraints_);

  // Get an array for the subset of l1 terms in the input vec.
  Eigen::Map<const Eigen::ArrayXd> l1_array(vec.data(), num_l1_residuals_);
  Eigen::Map<const Eigen::ArrayXd> inequality_array(
      vec.data() + num_l1_residuals_, num_inequality_constraints_);

  // Compute the L1 proximal operator on the L1 terms.
  output.head(num_l1_residuals_).array() =
      (l1_array - kappa).max(0.0) - (-l1_array - kappa).max(0.0);
  // Project the inequality constraints such that geq_mat * x - geq_vec > 0
  output.tail(num_inequality_constraints_).array() = inequality_array.max(0.0);
  return output;
}

}  // namespace theia
