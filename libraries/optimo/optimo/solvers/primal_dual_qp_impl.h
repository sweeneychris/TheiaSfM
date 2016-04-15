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

#ifndef OPTIMO_SOLVERS_PRIMAL_DUAL_QP_IMPL_H_
#define OPTIMO_SOLVERS_PRIMAL_DUAL_QP_IMPL_H_

#include <Eigen/Dense>

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::Dynamic;

// TODO(vfragoso): Allow to get the duals as well
template <typename Scalar, uint n, uint m, uint p>
TERMINATION_TYPE
PrimalDualQP<Scalar, n, m, p>::operator()(
    const typename PrimalDualQP<Scalar, n, m, p>::Params& params,
    Matrix<Scalar, n, 1>* x,
    Scalar* min_value) {
  if (!x || !min_value) return INVALID_ARGUMENTS;

  if (this->options.alpha_ < 0.0 ||
      this->options.alpha_ > 1.0 ||
      this->options.beta_ < 0.0 ||
      this->options.beta_ > 1.0 ||
      this->options.mu_ < 1.0) {
    return INVALID_ARGUMENTS;
  }

  // Allocate F_, and r_ and build the KKT system!
  const int l = n + m + p;  // KKT system dimensions
  F_ = Matrix<Scalar, Dynamic, Dynamic>::Zero(l, l);
  r_ = Matrix<Scalar, Dynamic, 1>::Zero(l);
  // Allocate our augmented vector of unkowns and duals
  Matrix<Scalar, Dynamic, 1> y(l);    // y = [x' lambdas' nus']'
  Matrix<Scalar, Dynamic, 1> ynt(l);  // The newton step
  Matrix<Scalar, Dynamic, 1> y_plus(l);  // y+ = y_ (for backtracking)
  Matrix<Scalar, Dynamic, 1> r_plus(l);  // r+ = r_ (for backtracking)
  // TODO(vfragoso): Check the initial state of the dual variables.
  // Leave the user to set the duals as it might not be feasible
  // or try to calculate a feasible lambda?
  y.setConstant(1.0);
  y.block(0, 0, n, 1) = *x;
  const auto& lambdas = y.block(n, 0, m, 1);
  //  const auto& nus = y.block(n + m, 0, p, 1);

  // Useful blocks
  Matrix<Scalar, Dynamic, 1> d(n);
  d.setConstant(0.0);
  d.block(0, 0, n, 1) = params.d;
  const auto& x_vec = y.block(0, 0, n, 1);
  // const auto& Q = F_.block(0, 0, n, n);
  const auto& Ain = F_.block(0, n, n, m).transpose();
  // const auto& Aeq = F_.block(0, n + m, n, p).transpose();

  // Check if initial point is feasible!
  Matrix<Scalar, Dynamic, 1> ineq_x = params.Ain*(*x) - params.bin;
  for (int i = 0; i < m; i++) {
    if (ineq_x(i) > 0.0) return INFEASIBLE_STARTING_POINT;
  }

  // TODO(vfragoso): Check if we can ignore equality constraints
  // Matrix<Scalar, Dynamic, 1> eq_x = params.Aeq*(*x) - params.beq;
  // for (int i = 0; i < p; i++) {
  //   if (fabs(eq_x(i)) > this->options.eps_feas_) {
  //     return INFEASIBLE_STARTING_POINT;
  //   }
  // }

  // Fill in F_
  // F=[Q_                 Ain'                  Aeq_';
  //    -diag(lambda)*Ain  -diag(fx)             zeros(m, size(Aeq,1));
  //    Aeq_               zeros(size(Aeq,1), m) zeros(size(Aeq,1), size(Aeq,1))
  // ]
  // fx are the inequalities evaluated: fx = Ain*x

  // First block row
  F_.block(0, 0, n, n) = params.Q;  // Placing Q
  F_.block(0, n, n, m) = params.Ain.transpose();  // Placing Ain
  F_.block(0, n + m, n, p) = params.Aeq.transpose();  // Placing Aeq
  // Third block row
  F_.block(n + m, 0, p, n) = params.Aeq;  // Placing Aeq

  for (int iter = 0; iter < this->options.max_iter_; iter++) {
    // 1. Build KKT System
    // Second block row
    F_.block(n, 0, m, n) = -1.0*(lambdas.asDiagonal())*Ain;
    Matrix<Scalar, m, 1> diag = Ain*x_vec - params.bin;
    F_.block(n, n, m, m) = -1.0*diag.asDiagonal();

    // Build residual vector
    // Dual residual
    // (Qx + d) + Ain.transpose()*lambdas + Aeq.transpose()*nus
    // Central residual
    // -diag(lambda)*(Ain*x - bin) - (1/t)*vec(ones)
    // Primal residual
    // Aeq*x - beq
    Scalar eta = -(diag.transpose()*lambdas)(0, 0);
    Scalar t = this->options.mu_*m/eta;
    const Scalar t_inv = 1.0/t;

    calculateResiduals(y, d, diag, params.beq, t_inv, &r_);
    // -- End of building kkt system

    // Should we stop?
    if (r_.block(0, 0, n, 1).norm() <= this->options.eps_feas_ &&
        r_.block(n + m, 0, p, 1).norm() <= this->options.eps_feas_ &&
        eta <= this->options.epsilon_) {
      x->block(0, 0, n, 1) = x_vec.block(0, 0, n, 1);
      // Eval QP obj
      *min_value =
          (0.5*x->transpose()*params.Q*(*x) + params.d.transpose()*(*x))(0, 0);
      return SOLVED;  // QP Solved!
    }

    // 2. Compute Newton step
    // ynt = F_.partialPivLu().solve(-r_);
    ynt = F_.householderQr().solve(-r_);

    // 3. Line search
    bool line_search_exit_flag;
    Scalar s = backTracking(t_inv, r_.norm(), ynt, y, d, params.beq, params.bin,
                            &y_plus, &r_plus, &line_search_exit_flag);

    // 4. Update
    y += s*ynt;
  }  // End of iteration loop

  return NOT_SOLVED;
}

template <typename Scalar, uint n, uint m, uint p>
inline void
PrimalDualQP<Scalar, n, m, p>::calculateResiduals(
    const Matrix<Scalar, Dynamic, 1>& y,
    const Matrix<Scalar, Dynamic, 1>& d,
    const Matrix<Scalar, Dynamic, 1>& diag,
    const Matrix<Scalar, p, 1>& beq,
    const Scalar t_inv,
    Matrix<Scalar, Dynamic, 1>* r) {
  const auto& Q = F_.block(0, 0, n, n);
  const auto& Ain = F_.block(0, n, n, m).transpose();
  const auto& Aeq = F_.block(0, n + m, n, p).transpose();

  const auto& x_vec = y.block(0, 0, n, 1);
  const auto& lambdas = y.block(n, 0, m, 1);
  const auto& nus = y.block(n + m, 0, p, 1);

  // Dual residual
  r->block(0, 0, n, 1) =
      Q*x_vec + d + Ain.transpose()*lambdas + Aeq.transpose()*nus;

  // Central residual
  r->block(n, 0, m, 1) = -t_inv*ones_ - lambdas.asDiagonal()*diag;

  // Primal residual
  r->block(n + m, 0, p, 1) = Aeq*x_vec - beq;
}

// TODO(vfragoso): Document me!!
// Algorithm from Boyd's convex optimization pg. 613
template <typename Scalar, uint n, uint m, uint p>
Scalar
PrimalDualQP<Scalar, n, m, p>::backTracking(
    const Scalar t_inv,
    const Scalar r_norm,
    const Matrix<Scalar, Dynamic, 1>& ynt,
    const Matrix<Scalar, Dynamic, 1>& y,
    const Matrix<Scalar, Dynamic, 1>& d,
    const Matrix<Scalar, p, 1>& beq,
    const Matrix<Scalar, m, 1>& bin,
    Matrix<Scalar, Dynamic, 1>* y_plus,
    Matrix<Scalar, Dynamic, 1>* r_plus,
    bool* exit_flag) {
  // Find s_max: Find from the newton step the largest entry in the vector.
  Scalar s = 1.00;  // Scalar to be returned
  const Scalar* ynt_data = ynt.data();
  for (uint i = 0; i < m; i++) {
    if (ynt_data[n + i] < static_cast<Scalar>(0.0)) {
      Scalar val = -y(n + i, 0) / ynt_data[n + i];
      if (val < s) s = val;
    }
  }

  // Regular backtracking
  s *= static_cast<Scalar>(0.99);
  Scalar rhs;
  Scalar lhs;
  uint ninf = 0;  // Number of infeasible inequalities.
  *y_plus = y;  // Copying y into y_plus

  // Useful matrix expressions
  const auto& Ain = F_.block(0, n, n, m).transpose();
  const auto& x_vec = y_plus->block(0, 0, n, 1);

  do {
    ninf = 0;
    *y_plus = y + s*ynt;
    Matrix<Scalar, m, 1> diag_fx = Ain*x_vec - bin;
    // Calculate residuals using y_plus
    calculateResiduals(*y_plus, d, diag_fx, beq, t_inv, r_plus);
    // Calculate lhs and rhs
    rhs = (static_cast<Scalar>(1.0) - this->options.alpha_*s)*r_norm;
    lhs = r_plus->norm();
    // Calculate number of infeasible inequalities
    for (int k = 0; k < m; k++) if (diag_fx(k) >= 0) ninf++;
    // Calculate new s
    s *= this->options.beta_;
  } while (lhs > rhs || ninf > 0);
  *exit_flag = true;
  s /= this->options.beta_;  // Since we always do an extra multiplication
  return s;
}
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_PRIMAL_DUAL_QP_IMPL_H_
