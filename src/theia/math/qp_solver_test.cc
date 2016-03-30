// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include <Eigen/Core>
#include <Eigen/LU>
#include "gtest/gtest.h"

#include "theia/math/qp_solver.h"
#include "theia/util/random.h"

namespace theia {

// A rigged QP minimization problem with a known output.
//
// 1/2 * x' * P * x + q' * x + r
//     [  5  -2  -1 ]
// P = [ -2   4   3 ]
//     [ -1   3   5 ]
//
// q = [ 2  -35  -47 ]^t
// r = 5
//
// Minimizing this unbounded problem should result in:
//   x = [ 3  5  7 ]^t
TEST(QPSolver, Unbounded) {
  static const double kTolerance = 1e-4;

  Eigen::MatrixXd P(3, 3);
  P << 5, -2, -1,
    -2, 4, 3,
    -1, 3, 5;
  Eigen::VectorXd q(3);
  q << 2, -35, -47;
  const double r = 5;

  QPSolver::Options options;
  options.max_num_iterations = 100;
  Eigen::SparseMatrix<double> P_sparse(P.sparseView());
  QPSolver qp_solver(options, P_sparse, q, r);
  Eigen::VectorXd solution;
  ASSERT_TRUE(qp_solver.Solve(&solution));

  // Verify the solution is near (3, 5, 7).
  const Eigen::Vector3d gt_solution(3, 5, 7);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(solution(i), gt_solution(i), kTolerance);
  }

  // Check that the residual is near optimal.
  const double residual =
      0.5 * solution.dot(P * solution) + solution.dot(q) + r;
  const double gt_residual =
      0.5 * gt_solution.dot(P * gt_solution) + gt_solution.dot(q) + r;
  EXPECT_NEAR(residual, gt_residual, kTolerance);
}

TEST(QPSolver, LooseBounds) {
  static const double kTolerance = 1e-4;

  Eigen::MatrixXd P(3, 3);
  P << 5, -2, -1,
    -2, 4, 3,
    -1, 3, 5;
  Eigen::VectorXd q(3);
  q << 2, -35, -47;
  const double r = 5;

  QPSolver::Options options;
  options.max_num_iterations = 100;
  Eigen::SparseMatrix<double> P_sparse(P.sparseView());
  QPSolver qp_solver(options, P_sparse, q, r);

  // Set a lower bound that should not affect the output.
  Eigen::VectorXd lower_bound(3);
  lower_bound << 0, 0, 0;
  qp_solver.SetLowerBound(lower_bound);

  // Set an upper bound that should not affect the output.
  Eigen::VectorXd upper_bound(3);
  upper_bound << 10, 10, 10;
  qp_solver.SetUpperBound(upper_bound);
  Eigen::VectorXd solution;

  ASSERT_TRUE(qp_solver.Solve(&solution));

  // Verify the solution is near (3, 5, 7).
  const Eigen::Vector3d gt_solution(3, 5, 7);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(solution(i), gt_solution(i), kTolerance);
  }

  // Check that the residual is near optimal.
  const double residual =
      0.5 * solution.dot(P * solution) + solution.dot(q) + r;
  const double gt_residual =
      0.5 * gt_solution.dot(P * gt_solution) + gt_solution.dot(q) + r;
  EXPECT_NEAR(residual, gt_residual, kTolerance);
}

TEST(QPSolver, TightBounds) {
  static const double kTolerance = 1e-4;

  Eigen::MatrixXd P(3, 3);
  P << 5, -2, -1,
    -2, 4, 3,
    -1, 3, 5;
  Eigen::VectorXd q(3);
  q << 2, -35, -47;
  const double r = 5;

  QPSolver::Options options;
  options.absolute_tolerance = 1e-8;
  options.relative_tolerance = 1e-8;
  Eigen::SparseMatrix<double> P_sparse(P.sparseView());
  QPSolver qp_solver(options, P_sparse, q, r);

  // Set a lower bound that constrains the output.
  Eigen::VectorXd lower_bound(3);
  lower_bound << 5, 7, 9;
  qp_solver.SetLowerBound(lower_bound);

  // Set an upper bound that constrains the output.
  Eigen::VectorXd upper_bound(3);
  upper_bound << 10, 12, 14;
  qp_solver.SetUpperBound(upper_bound);

  Eigen::VectorXd solution;
  ASSERT_TRUE(qp_solver.Solve(&solution));

  // Verify the solution is near (5, 7, 9).
  const Eigen::Vector3d gt_solution(5, 7, 9);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(solution(i), gt_solution(i), kTolerance);
  }

  // Check that the residual is near optimal.
  const double residual =
      0.5 * solution.dot(P * solution) + solution.dot(q) + r;
  const double gt_residual =
      0.5 * gt_solution.dot(P * gt_solution) + gt_solution.dot(q) + r;
  EXPECT_NEAR(residual, gt_residual, kTolerance);
}

TEST(QPSolver, InvalidBounds) {
  Eigen::MatrixXd P(3, 3);
  P << 5, -2, -1,
    -2, 4, 3,
    -1, 3, 5;
  Eigen::VectorXd q(3);
  q << 2, -35, -47;
  const double r = 5;

  QPSolver::Options options;
  Eigen::SparseMatrix<double> P_sparse(P.sparseView());
  QPSolver qp_solver(options, P_sparse, q, r);

  // Set the upper bound as the lower bound.
  Eigen::VectorXd lower_bound(3);
  lower_bound << 5, 7, 9;
  qp_solver.SetUpperBound(lower_bound);

  // Set the lower bound as the upper bound, making the valid solution space
  // non-existant.
  Eigen::VectorXd upper_bound(3);
  upper_bound << 10, 12, 14;
  qp_solver.SetLowerBound(upper_bound);

  Eigen::VectorXd solution;
  EXPECT_FALSE(qp_solver.Solve(&solution));
}

}  // namespace theia
