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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#ifndef THEIA_MATH_CONSTRAINED_L1_SOLVER_H_
#define THEIA_MATH_CONSTRAINED_L1_SOLVER_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <glog/logging.h>

#include "theia/math/matrix/sparse_cholesky_llt.h"

namespace theia {

// This class allows for solving linear L1 minimization problems under
// inequality constraints. These problems are of the form:
//
//   minimize    ||Ax - b||_1
//   subject to:   Cx - d > 0
//
// This convex problem is solved with the alternating direction method of
// multipliers (ADMM). We do not implement the "less than" constraint because
// the functions that utilize this solver do not need them.
//
// ADMM can be much faster than interior point methods but convergence may be
// slower. Generally speaking, ADMM solvers converge to good solutions in only a
// few number of iterations, but can spend many iterations subsuquently refining
// the solution to optain the global optimum. The speed improvements are because
// the linear system only needs to be factorized (e.g., by Cholesky
// decomposition) once, as opposed to every iteration.
class ConstrainedL1Solver {
 public:
  struct Options {
    int max_num_iterations = 1000;
    // Rho is the augmented Lagrangian parameter for the L1 minimization.
    double rho = 10.0;
    // Alpha is the over-relaxation parameters (typically between 1.0 and 1.8).
    double alpha = 1.2;

    // Stopping criteria.
    double absolute_tolerance = 1e-4;
    double relative_tolerance = 1e-2;
  };

  // The linear system along with the equality and inequality constraints.
  ConstrainedL1Solver(const Options& options,
                      const Eigen::SparseMatrix<double>& A,
                      const Eigen::VectorXd& b,
                      const Eigen::SparseMatrix<double>& geq_mat,
                      const Eigen::VectorXd& geq_vec);

  // Solve the constrained L1 minimization above.
  void Solve(Eigen::VectorXd* solution);

 private:
  // This method is used for the z-update, which is conveniently an element-wise
  // update. For the terms in vec corresponding to the L1 minimization, we
  // update
  // the values with the L1 proximal mapping (Shrinkage) operator. The terms
  // corresponding to the inequality constraints are constrained to be greater
  // than zero as vec = max(vec, 0).
  Eigen::VectorXd ModifiedShrinkage(const Eigen::VectorXd& vec,
                                    const double kappa);

  const Options options_;
  const int num_l1_residuals_;
  const int num_inequality_constraints_;

  // Matrix A where || Ax - b ||_1 is the problem we are solving.
  Eigen::SparseMatrix<double> A_;
  Eigen::VectorXd b_;

  // Cholesky linear solver. Since our linear system will be a SPD matrix we can
  // utilize the Cholesky factorization.
  SparseCholeskyLLt linear_solver_;
};

}  // namespace theia

#endif  // THEIA_MATH_CONSTRAINED_L1_SOLVER_H_
