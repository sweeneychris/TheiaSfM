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

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <glog/logging.h>

namespace theia {

// These are template overrides that allow the sparse linear solvers to work
// with sparse or dense matrices. The sparseView() method is not implemented for
// Eigen::SparseMatrix.
namespace l1_solver_internal {

inline void AnalyzePattern(
    const Eigen::SparseMatrix<double>& spd_mat,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >* linear_solver) {
  linear_solver->analyzePattern(spd_mat);
}

inline void AnalyzePattern(
    const Eigen::MatrixXd& spd_mat,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >* linear_solver) {
  linear_solver->analyzePattern(spd_mat.sparseView());
}

inline void Factorize(
    const Eigen::SparseMatrix<double>& spd_mat,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >* linear_solver) {
  linear_solver->factorize(spd_mat);
}

inline void Factorize(
    const Eigen::MatrixXd& spd_mat,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >* linear_solver) {
  linear_solver->factorize(spd_mat.sparseView());
}

}  // namespace l1_solver_internal

// An L1 norm approximation solver. This class will attempt to solve the
// problem: || A * x - b || under L1-norm (as opposed to L2 i.e. "least-squares"
// norm). This problem can be cast as a simple linear program which, in turn, is
// actually just a simple set of weighted least squares problems. We use a
// MatrixType template type so that dense or sparse matrices can be used,
// however, they must be a double type.
//
// The solution strategy comes from the book "Convex Optimization" by Boyd and
// Vandenberghe: http://web.stanford.edu/~boyd/cvxbook/
// Section 11.8.2 gives an implementation of how to solve the L1 norm
// approximation problem with Newton steps, which is the method used here.
template <class MatrixType>
class L1Solver {
 public:
  struct Options {
    // If the primal dual gap is smaller than this then the algorithm stops.
    double duality_gap_tolerance = 1e-3;
    int max_num_iterations = 100;

    // Parameters for convergence and backtracking.
    double alpha = 0.01;
    double beta = 0.5;
    double mu = 20;
  };

  L1Solver(const Options& options, const MatrixType& mat)
      : options_(options), a_(mat) {
    // Analyze the sparsity pattern once. Only the values of the entries will be
    // changed with each iteration.
    const MatrixType spd_mat = a_.transpose() * a_;
    l1_solver_internal::AnalyzePattern(spd_mat, &linear_solver_);
    CHECK_EQ(linear_solver_.info(), Eigen::Success);
  }

  // Solves ||Ax - b||_1 for the optimial L1 solution given an initial guess for
  // x. To solve this we introduce an auxillary variable y such that the
  // solution to:
  //        min   1 * y
  //   s.t. [  A   -I ] [ x ] < [  b ]
  //        [ -A   -I ] [ y ]   [ -b ]
  // which is an equivalent linear program.
  void Solve(const Eigen::VectorXd& rhs, Eigen::VectorXd* solution) {
    x_ = solution;
    rhs_ = rhs;

    ax_ = a_ * (*x_);
    const Eigen::VectorXd initial_l2_residual = ax_ - rhs_;
    y_ = 0.95 * (initial_l2_residual.array().abs()) +
         0.1 * (initial_l2_residual.array().abs()).maxCoeff();

    // Initialize the primal_penalty and dual variables.
    primal_penalty1_ = initial_l2_residual - y_;
    primal_penalty2_ = -initial_l2_residual - y_;
    lambda1_ = -primal_penalty1_.array().inverse();
    lambda2_ = -primal_penalty2_.array().inverse();
    atv_ = a_.transpose() * (lambda1_ - lambda2_);

    for (int i = 0; i < options_.max_num_iterations; i++) {
      const double surrogate_duality_gap =
          -(primal_penalty1_.dot(lambda1_) + primal_penalty2_.dot(lambda2_));

      // TODO(cmsweeney): Check the dual residual for convergence.
      if (surrogate_duality_gap <= options_.duality_gap_tolerance) {
        VLOG(1) << "Converged in " << i + 1 << " iterations.";
        break;
      }
      const double tau = options_.mu * rhs_.size() / surrogate_duality_gap;

      // Solve for the direction of the newton step. For L1 minimization this is
      // a special-case which is more simple than the general LP.
      if (!ComputeNewtonStep(tau)) {
        LOG(WARNING) << "Could not compute Newton step.";
        return;
      }

      // Compute the maximum step size to remain a feasible solution.
      ComputeStepSize(tau);
    }
  }

 private:
  // Determines the primal-dual search direction from the linear program.
  bool ComputeNewtonStep(const double tau) {
    const double inv_tau = 1.0 / tau;
    const Eigen::VectorXd x_quotient = lambda1_.cwiseQuotient(primal_penalty1_);
    const Eigen::VectorXd y_quotient = lambda2_.cwiseQuotient(primal_penalty2_);
    const Eigen::VectorXd sig1 = -x_quotient - y_quotient;
    const Eigen::VectorXd sig2 = x_quotient - y_quotient;
    const Eigen::VectorXd sigx = sig1 - sig2.cwiseAbs2().cwiseQuotient(sig1);
    const Eigen::VectorXd w1 =
        -inv_tau *
        (a_.transpose() * (-primal_penalty1_.array().inverse() +
                             primal_penalty2_.array().inverse()).matrix());
    const Eigen::VectorXd w2 =
        -1.0 - inv_tau * ((primal_penalty1_.array().inverse() +
                           primal_penalty2_.array().inverse())).array();

    const Eigen::VectorXd w1p =
        w1 - a_.transpose() * ((sig2.cwiseQuotient(sig1)).cwiseProduct(w2));
    const MatrixType lhs = a_.transpose() * sigx.asDiagonal() * a_;

    // Factorize the matrix based on the current linear system. If factorization
    // fails, return false.
    l1_solver_internal::Factorize(lhs, &linear_solver_);
    if (linear_solver_.info() != Eigen::Success) {
      LOG(WARNING) << "Failed to compute a Sparse Cholesky factorization for "
                      "Simplicial LDLT.";
      return false;
    }

    dx_ = linear_solver_.solve(w1p);
    CHECK_EQ(linear_solver_.info(), Eigen::Success);

    adx_ = a_ * dx_;

    dy_ = (w2 - sig2.cwiseProduct(adx_)).cwiseQuotient(sig1);
    dlambda1_ =
        -(lambda1_.cwiseQuotient(primal_penalty1_)).cwiseProduct(adx_ - dy_) -
        lambda1_ - (inv_tau * primal_penalty1_.array().inverse()).matrix();
    dlambda2_ =
        (lambda2_.cwiseQuotient(primal_penalty2_)).cwiseProduct(adx_ + dy_) -
        lambda2_ - (inv_tau * primal_penalty2_.array().inverse()).matrix();
    return true;
  }

  // Computes a step size to ensure that the solution will lie within the
  // feasible region, i.e., lambda1, lambda2 > 0 and primal_penalty1,
  // primal_penalty2 < 0.
  double ComputeFeasibleStep() {
    double step_size = 1.0;

    // Ensure lambda1 and lambda 2 > 0.
    for (int i = 0; i < dlambda1_.size(); i++) {
      if (dlambda1_(i) < 0) {
        const double new_step_size = -lambda1_(i) / dlambda1_(i);
        if (new_step_size < step_size) {
          step_size = new_step_size;
        }
      }
      if (dlambda2_(i) < 0) {
        const double new_step_size = -lambda2_(i) / dlambda2_(i);
        if (new_step_size < step_size) {
          step_size = new_step_size;
        }
      }
    }

    // Ensure that primal_penalty1, primal_penalty2 < 0.
    for (int i = 0; i < primal_penalty1_.size(); i++) {
      if (adx_(i) - dy_(i) > 0) {
        const double new_step_size = -primal_penalty1_(i) / (adx_(i) - dy_(i));
        if (new_step_size < step_size) {
          step_size = new_step_size;
        }
      }
      if (-adx_(i) - dy_(i) > 0) {
        const double new_step_size = -primal_penalty2_(i) / (-adx_(i) - dy_(i));
        if (new_step_size < step_size) {
          step_size = new_step_size;
        }
      }
    }
    step_size *= 0.99;
    return step_size;
  }

  // Computes the step size using backtracking given the direction of the newton
  // step.
  void ComputeStepSize(const double tau) {
    static const int kMaxBacktrackIterations = 32;
    const double inv_tau = 1.0 / tau;

    const Eigen::VectorXd atdv = a_.transpose() * (dlambda1_ - dlambda2_);
    double rdual =
        ((-lambda1_ - lambda2_).array() + 1.0).matrix().squaredNorm();

    Eigen::VectorXd temp1 =
        (-lambda1_.cwiseProduct(primal_penalty1_)).array() + inv_tau;
    Eigen::VectorXd temp2 =
        (-lambda2_.cwiseProduct(primal_penalty2_)).array() + inv_tau;
    const double norm_residual = std::sqrt(
        atv_.squaredNorm() + rdual + temp1.squaredNorm() + temp2.squaredNorm());

    // Make sure that the step size is feasible i.e. lambda1, lambda2 > 0 and
    // primal_penalty1, primal_penalty2 > 0.
    double step_size = ComputeFeasibleStep();

    Eigen::VectorXd x_p, y_p, lambda1_p, lambda2_p, axp, atvp;
    for (int i = 0; i < kMaxBacktrackIterations; i++) {
      x_p = *x_ + step_size * dx_;
      y_p = y_ + step_size * dy_;

      axp = ax_ + step_size * adx_;
      atvp = atv_ + step_size * atdv;

      lambda1_p = lambda1_ + step_size * dlambda1_;
      lambda2_p = lambda2_ + step_size * dlambda2_;

      primal_penalty1_ = axp - rhs_ - y_p;
      primal_penalty2_ = -axp + rhs_ - y_p;

      rdual = ((-lambda1_p - lambda2_p).array() + 1.0).matrix().squaredNorm();
      temp1 = (-lambda1_.cwiseProduct(primal_penalty1_)).array() + inv_tau;
      temp2 = (-lambda2_.cwiseProduct(primal_penalty2_)).array() + inv_tau;
      const double step_norm =
          std::sqrt(atvp.squaredNorm() + rdual + temp1.squaredNorm() +
                    temp2.squaredNorm());
      step_size *= options_.beta;
      if (step_norm <= (1.0 - options_.alpha * step_size) * norm_residual) {
        break;
      }
    }

    *x_ = x_p;
    y_ = y_p;
    lambda1_ = lambda1_p;
    lambda2_ = lambda2_p;
    ax_ = axp;
    atv_ = atvp;
  }

  Options options_;

  // Solution vector.
  Eigen::VectorXd* x_;

  // Matrix A where || Ax - b ||_1 is the problem we are solving.
  MatrixType a_;

  // rhs corresponds to the vector b, and y is the auxillary variable that we
  // minimize for the linear program:
  //   minimize y s.t.
  //    Ax - b < y
  //   -Ax + b < y
  Eigen::VectorXd y_, rhs_;

  // Primal penalties correspond to the inverse of lambda, the dual variables.
  Eigen::VectorXd primal_penalty1_, primal_penalty2_;
  Eigen::VectorXd lambda1_, lambda2_;

  // Derivitives computed by the Newton step.
  Eigen::VectorXd dx_, dy_, dlambda1_, dlambda2_;

  // Pre-computed values A * x, A * dx, A^t * (lambda1 - lambda2).
  Eigen::VectorXd ax_, adx_, atv_;

  // Cholesky linear solver. Since our linear system will be a SPD matrix we can
  // utilize the Cholesky factorization.
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > linear_solver_;
};

}  // namespace theia

#endif  // THEIA_MATH_L1_SOLVER_H_
