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

#ifndef THEIA_MATH_L1_SOLVER_H_
#define THEIA_MATH_L1_SOLVER_H_

#include <Eigen/Core>
#include <Eigen/SPQRSupport>

namespace theia {

// An L1 norm approximation solver. This class will attempt to solve the
// problem: || A * x - b || under L1-norm (as opposed to L2 i.e. "least-squares"
// norm). This problem can be cast as a simple linear program which, in turn, is
// actually just a simple set of weighted least squares problems. We use a
// MatrixBase template type so that dense or sparse matrices can be used,
// however, they must be a double type.
//
// The solution strategy comes from the book "Convex Optimization" by Boyd and
// Vandenberghe: http://web.stanford.edu/~boyd/cvxbook/
// Section 11.8.2 gives an implementation of how to solve the L1 norm
// approximation problem with Newton steps, which is the method used here.
template <class MatrixBase>
class L1Solver {
  struct Options {
    // The minimum accuracy of the solution.
    double epsilon = 1e-8;
    double mu = 10.0;
    double alpha = 0.01;
    double beta = 0.5;
    int max_num_iterations = 100;
  };

  L1Solver(const Options& options, const MatrixBase& mat)
      : options_(options), mat_(mat) {
    spqr_.compute(mat_);
    CHECK_EQ(spqr_.info(), Eigen::Success);

    // TODO(cmsweeney): how to initialize t_?
  }

  // Solves ||Ax - b||_1 for the optimial L1 solution given an initial guess for
  // x.
  void Solve(const Eigen::VectorXd& initial_x,
             const Eigen::VectorXd& rhs,
             Eigen::VectorXd* solution) {
    Eigen::VectorXd y;
    // TODO(cmsweeney): How to initialize y?
    // TODO(cmsweeney): How to initialize t_?

    for (int i = 0; i < options_.max_num_iterations; i++) {
      ComputeNewtonStep(*solution, y, rhs);
      const double step_size = ComputeStepSize(*solution);
      *solution += step_size * dx_;
      y += step_size * dy_;

      // Determing if the step size is small enough to indicate convergence.
      if (solution->size() * 2 / t_ < options_.epsilon) {
        VLOG(1) << "Converged in " << i + 1 << " iterations.";
        break;
      }

      t_ *= options_.mu;
    }
  }

 private:
  // From 11.8.2 in the "Convex Optimization" book.
  void ComputeNewtonStep(const Eigen::VectorXd& x,
                         const Eigen::VectorXd& y,
                         const Eigen::VectorXd& b) {
    // D1 = diag(b - Ax + y)^-2
    const Eigen::VectorXd inv_positive_solution =
        (b - mat_ * x + y).array().inverse();
    const Eigen::VectorXd d1 = inv_positive_solution.array().square();

    // D2 = diag(-b + Ax + y)^-2
    const Eigen::VectorXd inv_negative_solution =
        (-b - mat_ * x + y).array().inverse();
    const Eigen::VectorXd d2 = inv_negative_solution.array().square();

    // D = 4 * D1 * D2 * (D1 + D2)^-1 = 2 * (diag(y)^2 + diag(b - Ax)^2)^-1
    const Eigen::VectorXd d = 4.0 * d1.asDiagonal() * d2.asDiagonal() *
                              ((d1 + d2).asDiagonal().inverse());

    // g1 = diag(b - Ax + y)^-1 - diag(-b + Ax + y)^-1
    const Eigen::VectorXd g1 = inv_positive_solution - inv_negative_solution;

    // g2 = t - diag(b - Ax + y)^-1 - diag(-b + Ax + y)^-1
    const Eigen::VectorXd g2 =
        t_ - inv_positive_solution.array() - inv_negative_solution.array();

    // g = g1 + (D1 - D2) * (D1 + D2)^-1 * g2
    const Eigen::VectorXd g =
        g1 + (d1 - d2).asDiagonal() * ((d1 + d2).asDiagonal().inverse()) * g2;

    // Solve ||A * dx - D^-1 * g||_2.
    dx_ = spqr_.solve(d.asDiagonal().inverse() * g);
    CHECK_EQ(spqr_.info(), Eigen::Success);

    // dy = (D1 + D2)^-1 * (-g2 + (D1 - D2) * A * dx).
    dy_ = (d1 + d2).asDiagonal().inverse() *
         (-g2 + (d1 - d2).asDiagonal() * mat_ * dx);
  }

  // Computes the step size using backtracking given the direction of the newton
  // step.
  double ComputeStepSize(const Eigen::VectorXd& x) {
    // TODO(cmsweeney): How to intialize this?
    double step_size = 0.99;

    const double fx = (mat_ * x).squaredNorm();
    double fx_plus_dx = (mat * (x + step_size * dx)).norm();
    while (fx_plus_dx <= (1.0 - step_size * options_.alpha) * fx) {
      step_size *= options_.beta;
      fx_plus_dx = (mat * (x + step_size * dx)).norm();
    }
    return step_size;
  }

  Options options_;
  Eigen::SPQR<MatrixBase> spqr_;
  Eigen::VectorXd dx_, dy_;

  // TODO(cmsweeney): How to initialize this??
  double t_ = 0.0;
  const int max_iterations_ = 100;
};

}  // namespace theia

#endif  // THEIA_MATH_L1_SOLVER_H_
