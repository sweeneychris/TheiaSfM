// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#ifndef THEIA_MATH_QP_SOLVER_H_
#define THEIA_MATH_QP_SOLVER_H_

#include <Eigen/Core>

#include "theia/math/matrix/sparse_cholesky_llt.h"

namespace theia {

// An Quadratic Program (QP) solver. This class will attempt to
//    minimize: 1/2 * x' * P * x + q' * x + r
//    subject to: lb <= x <= ub
//
// This problem can be solved with the alternating direction method of
// multipliers (ADMM) as a least unsquared deviations minimizer. A full
// description of the method, including how to use ADMM for L1 minimization can
// be found in "Distributed Optimization and Statistical Learning via the
// Alternating Direction Method of Multipliers" by Boyd et al, Foundations and
// Trends in Machine Learning (2012). The paper can be found at:
//   https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
//
// ADMM can be much faster than interior point methods but convergence may be
// slower. Generally speaking, ADMM solvers converge to good solutions in only a
// few number of iterations, but can spend many iterations subsuquently refining
// the solution to optain the global optimum. The speed improvements are because
// the matrix A only needs to be factorized (by Cholesky decomposition) once, as
// opposed to every iteration.
//
// This implementation is based off of the code found at:
//   https://web.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
class QPSolver {
 public:
  struct Options {
    int max_num_iterations = 1000;
    // Rho is the augmented Lagrangian parameter.
    double rho = 1.0;
    // Alpha is the over-relaxation parameter (typically between 1.0 and 1.8).
    double alpha = 1.0;

    double absolute_tolerance = 1e-6;
    double relative_tolerance = 1e-4;
  };

  // Set Q, p, and r according to the notation above.
  QPSolver(const Options& options,
           const Eigen::SparseMatrix<double>& P,
           const Eigen::VectorXd& q,
           const double r);

  // Set the maximum number of iterations for the solver.
  void SetMaxIterations(const int max_iterations);

  // Set upper and lower bounds for the solution.
  void SetUpperBound(const Eigen::VectorXd& ub);
  void SetLowerBound(const Eigen::VectorXd& lb);

  // Solve the quadratic program.
  bool Solve(Eigen::VectorXd* solution);

 private:
  Options options_;

  // Matrix P, q, double r where || 1/2 * x * P * x + q' * x + b ||_2 is the
  // problem we are solving.
  const Eigen::SparseMatrix<double>& P_;
  const Eigen::VectorXd& q_;
  const double r_;

  // Lower and upper bounds.
  Eigen::ArrayXd lb_, ub_;

  // Cholesky linear solver. Since our linear system will be a SPD matrix we can
  // utilize the Cholesky factorization.
  SparseCholeskyLLt linear_solver_;
};

}  // namespace theia

#endif  // THEIA_MATH_QP_SOLVER_H_
