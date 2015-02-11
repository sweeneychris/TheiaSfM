// Copyright (C) 2014  Victor Fragoso <vfragoso@cs.ucsb.edu>
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
//     * Neither the name of the University of California, Santa Barbara nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL VICTOR FRAGOSO BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef OPTIMO_SOLVERS_PRIMAL_DUAL_LP_IMPL_H_
#define OPTIMO_SOLVERS_PRIMAL_DUAL_LP_IMPL_H_

#include <map>

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::DiagonalMatrix;

// Dense Primal Dual Implementation
template <typename Scalar, uint n, uint m, uint p>
TERMINATION_TYPE
PrimalDualLP<Scalar, n, m, p>::operator()(const Params& params,
                                          Matrix<Scalar, Dynamic, 1>* y,
                                          Scalar* min_value) {
  const int l = n + m + p;  // Size of F_
  // Check arguments
  if (!y || !min_value || y->rows() != l) return INVALID_ARGUMENTS;

  // Checking solver options
  if (this->options.alpha_ < 0.0 ||
      this->options.alpha_ > 1.0 ||
      this->options.beta_ < 0.0 ||
      this->options.beta_ > 1.0 ||
      this->options.mu_ < 1.0) {
    return INVALID_ARGUMENTS;
  }

  F_ = Matrix<Scalar, Dynamic, Dynamic>::Zero(l, l);
  y->block(n, 0, m + p, 1).setConstant(1.0);  // Setting duals to 1
  r_ = Matrix<Scalar, Dynamic, 1>(l);
  Matrix<Scalar, Dynamic, 1> ynt(l);  // Newton Step
  y_plus_ = Matrix<Scalar, Dynamic, 1>(l);
  r_plus_ = Matrix<Scalar, Dynamic, 1>(l);

  Scalar eta;  // Surrogate gap
  Scalar t;  // Term associated to log-barrier
  Scalar s;  // Step size

  // Checking if F_ really had data in it
  if (params.Aeq.sum() == 0.0 || params.Ain.sum() == 0.0) {
    return INVALID_ARGUMENTS;
  }

  F_.block(0, n, n, m) = params.Ain.transpose();  // Placing Ain
  F_.block(0, n + m, n, p) = params.Aeq.transpose();  // Placing Aeq
  F_.block(n + m, 0, p, n) = params.Aeq;  // Placing Aeq

  // Check if initial point is feasible!
  const auto& x = y->block(0, 0, n, 1);
  Matrix<Scalar, Dynamic, 1> ineq_x = params.Ain*x - params.bin;
  for (int i = 0; i < m; i++) {
    if (ineq_x(i) > 0.0) return INFEASIBLE_STARTING_POINT;
  }

  // TODO(vfragoso): Do we really need to check Equality?
  // Matrix<Scalar, Dynamic, 1> eq_x = params.Aeq*x - params.beq;
  // for (int i = 0; i < p; i++) {
  //   if (fabs(eq_x(i)) > this->options.eps_feas_) {
  //     return INFEASIBLE_STARTING_POINT;
  //   }
  // }

  // Useful blocks
  const auto& lambdas = y->block(n, 0, m, 1);
  const auto& nus = y->block(n + m, 0, p, 1);

  // Allocating space for residuals
  double rp_norm, rd_norm;
  for (int iter = 0; iter < this->options.max_iter_; iter++) {  // Main loop
    // 1. Determine t and build/update KKT system
    buildKKTSystem(params, *y, this->options.mu_, &eta, &t, &r_,
                   &rd_norm, &rp_norm);

    // Should we stop?
    if (rd_norm <= this->options.eps_feas_ &&
        rp_norm <= this->options.eps_feas_ &&
        eta <= this->options.epsilon_) {
      *min_value = params.c.dot(x);
      return SOLVED;  // LP Solved!
    }
    // 2. Compute Newton step
    ynt = F_.householderQr().solve(-r_);

    // 3. Line Search
    s = backTracking(static_cast<Scalar>(1.0/t), r_.norm(),
                     ynt, *y, params);

    // 4. Update
    *y += s*ynt;
  }

  return MAX_ITERATIONS;
}

template <typename Scalar, uint n, uint m, uint p>
void
PrimalDualLP<Scalar, n, m, p>::buildKKTSystem(
    const Params& params,
    const Matrix<Scalar, Dynamic, 1>& y,
    const Scalar mu,
    Scalar* eta,
    Scalar* t,
    Matrix<Scalar, Dynamic, 1>* residuals,
    Scalar* rd_norm,
    Scalar* rp_norm) {
  // Build Second Block-Row of F
  // F=[zeros(n, n)       Ain'                  Aeq';
  //    -diag(lambda)*Ain -diag(fx)             zeros(m, size(Aeq,1));
  //    Aeq               zeros(size(Aeq,1), m) zeros(size(Aeq,1), size(Aeq,1))]
  // This implementation expects that Ain and Aeq are stored in the first
  // block-row. The first and third rows are automatically set when the user
  // sets Aeq and Ain. The second row is calculated below.

  const auto& x = y.block(0, 0, n, 1);  // Primal variables
  const auto& lambdas = y.block(n, 0, m, 1);  // Dual variables (ineq)
  const auto& nu = y.block(n + m, 0, p, 1);  // Dual variables (eq)
  Matrix<Scalar, Dynamic, 1> fx = params.Ain*x - params.bin;
  DiagonalMatrix<Scalar, Dynamic, Dynamic> diag_lambdas(lambdas);
  F_.block(n, 0, m, n) =
      static_cast<Scalar>(-1.0)*diag_lambdas*params.Ain;  // first block
  F_.block(n, n, m, m) =
      static_cast<Scalar>(-1.0)*fx.asDiagonal();
  
  *eta = static_cast<Scalar>(-fx.dot(lambdas));
  *t = mu * m / *eta;
  Scalar t_inv = static_cast<Scalar>(1.0 / *t);
  calculateResiduals(params, y, diag_lambdas, fx, t_inv, residuals);
  *rd_norm = residuals->block(0, 0, n, 1).norm();
  *rp_norm = residuals->block(n + m, 0, p, 1).norm();
}

template <typename Scalar, uint n, uint m, uint p>
void
PrimalDualLP<Scalar, n, m, p>::calculateResiduals(
    const Params& params,
    const Matrix<Scalar, Dynamic, 1>& y,
    const DiagonalMatrix<Scalar, Dynamic, Dynamic>& diag,
    const Matrix<Scalar, Dynamic, 1>& fx,
    const Scalar t_inv,
    Matrix<Scalar, Dynamic, 1>* residuals) {
  // Useful blocks
  const auto& x = y.block(0, 0, n, 1);  // Primal variables
  const auto& lambdas = y.block(n, 0, m, 1);  // Dual variables (ineq)
  const auto& nu = y.block(n + m, 0, p, 1);  // Dual variables (eq)

  // Build Residuals
  // Dual residuals
  const int l = n + m + p;
  residuals->block(0, 0, n, 1) =
      params.Ain.transpose()*lambdas + params.Aeq.transpose()*nu + params.c;
  // Central residuals
  residuals->block(n, 0, m, 1) =
       -t_inv*Matrix<Scalar, Dynamic, 1>::Ones(m) - diag*fx;
  // Primal residuals
  residuals->block(n + m, 0, p, 1) = params.Aeq*x - params.beq;
}

template <typename Scalar, uint n, uint m, uint p>
Scalar PrimalDualLP<Scalar, n, m, p>::backTracking(
    const double t_inv,
    const double r_norm,
    const Matrix<Scalar, Dynamic, 1>& ynt,
    const Matrix<Scalar, Dynamic, 1>& y,
    const Params& params) {
  // Find s_max
  Scalar s = 1.00;
  const Scalar* ynt_data = ynt.data();
  const Scalar* y_data = y.data();
  for (uint i = 0; i < m; i++) {
    if (ynt_data[n + i] < static_cast<Scalar>(0.0)) {
      Scalar val = -y_data[n + i]/ynt_data[n + i];
      if (val < s) s = val;
    }
  }

  // Do regular backtracking
  const int l = n + m + p;
  const Scalar zero = static_cast<Scalar>(0.0);
  s *= static_cast<Scalar>(0.99);
  double rhs;
  double lhs;
  uint ninf = 0;  // Number of infeasible ineq.

  // Useful blocks
  const auto& x_vec = y_plus_.block(0, 0, n, 1);  // x plus
  const auto& lambdas_vec = y_plus_.block(n, 0, m, 1);  // Lambdas plus

  do {
    ninf = 0;
    y_plus_ = y + s*ynt;
    // Calculate Residuals w/ y_plus
    Matrix<Scalar, Dynamic, 1> fx = params.Ain*x_vec - params.bin;
    DiagonalMatrix<Scalar, Dynamic, Dynamic> diag_lambdas(lambdas_vec);
    calculateResiduals(params, y_plus_, diag_lambdas, fx, t_inv, &r_plus_);
    // Calculate rhs and lhs
    rhs = static_cast<Scalar>(1.0 - this->options.alpha_*s)*r_norm;
    lhs = r_plus_.norm();
    // Calculate number of infeasible inequalities
    for (int k = 0; k < m; k++) if (fx(k) >= 0) ninf++;
    s *= this->options.beta_;
  } while (lhs > rhs || ninf > 0);
  s /= this->options.beta_;  // Since we always do an extra multiplication
  return s;
}
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_PRIMAL_DUAL_LP_IMPL_H_
