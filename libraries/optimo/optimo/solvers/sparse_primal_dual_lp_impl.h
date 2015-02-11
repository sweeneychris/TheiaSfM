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

#ifndef OPTIMO_SOLVERS_SPARSE_PRIMAL_DUAL_LP_IMPL_H_
#define OPTIMO_SOLVERS_SPARSE_PRIMAL_DUAL_LP_IMPL_H_

#include <unordered_map>

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::SparseMatrix;
using Eigen::Dynamic;
using Eigen::DiagonalMatrix;

// Implementation of LP solvers
template <typename Scalar>
TERMINATION_TYPE
SparsePrimalDualLP<Scalar>::operator()(Params* params,
                                       Matrix<Scalar, Dynamic, 1>* y,
                                       Scalar* min_value) {
    // Check for valid Primal Dual Newton parameters
  if (this->options.alpha_ < 0.0 || this->options.alpha_ > 1.0 ||
      this->options.beta_ < 0.0 || this->options.beta_ > 1.0 ||
      this->options.mu_ < 1.0) {
    return INVALID_ARGUMENTS;
  }

  // Check input
  if (!params || !y || !min_value || y->rows() != params->l_) {
    return INVALID_ARGUMENTS;
  }

  // Set pointers to get easier access and dimensions
  y_ = y;
  F_ = &params->F_;
  beq_ = &params->beq;
  bin_ = &params->bin;
  c_ = &params->c;
  n_ = params->n_;
  m_ = params->m_;
  p_ = params->p_;
  l_ = params->l_;

  // Check that F was filled
  if (!F_ || !beq_ || !bin_ || !c_) return INVALID_ARGUMENTS;

  if (isF_empty()) {
    return INVALID_ARGUMENTS;
  }

  int iter = 0;
  Matrix<Scalar, Eigen::Dynamic, 1> ynt(l_);  // Newton Step
  y_->block(n_, 0, m_ + p_, 1).setConstant(static_cast<Scalar>(1.0));
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y_plus(l_);

  // Testing Feasibility of initial point
  if (!feasibleStartingPoint(*y)) {
    return INFEASIBLE_STARTING_POINT;
  }

  Scalar eta;  // Surrogate gap
  Scalar t;  // "log barrier" term for central residual
  Scalar s;  // step size (linesearch)

  // Allocating space for residuals
  Matrix<Scalar, Dynamic, 1> residuals(l_, 1);
  Matrix<Scalar, Dynamic, 1> residuals_plus(l_, 1);
  Scalar rp_norm, rd_norm;

  Eigen::SPQR<SparseMatrix<Scalar> > lu_solver;
  while (iter++ < this->options.max_iter_) {  // Main Loop
    // 1. Determine t and build/update KKT system
    buildKKTSystem(this->options.mu_, &eta, &t, &residuals, &rd_norm, &rp_norm);

    if (!F_->isCompressed()) F_->makeCompressed();

    // Should we stop?
    if (rd_norm <= this->options.eps_feas_ &&
        rp_norm <= this->options.eps_feas_ &&
        eta <= this->options.epsilon_) {
      *min_value = params->c.dot(y->block(0, 0, n_, 1));
      return SOLVED;  // LP Solved!
    }

    // 2. Compute Newton step
    // Solve KKT system
    lu_solver.compute(*F_);
    if (lu_solver.info() != Eigen::Success) {
      return mapResult(lu_solver);
    }

    // Compute Newton Step
    ynt = -1.0*lu_solver.solve(residuals);
    if (lu_solver.info() != Eigen::Success) {
      return mapResult(lu_solver);
    }

    // 3. Line Search
    if (!backTracking(1.0/t, residuals.norm(), ynt,
                      &y_plus, &residuals_plus, &s)) {
      return NO_CONVERGENCE;
    }

    // 4. Update
    *y_ += s*ynt;
  }

  return SOLVED;
}

template <typename Scalar>
void SparsePrimalDualLP<Scalar>::buildKKTSystem(
    const double mu,
    double* eta,
    double* t,
    Matrix<Scalar, Dynamic, 1>* residuals,
    double* rd_norm,
    double* rp_norm) {
  // Build Second Block-Row of F
  // F=[zeros(n, n)       Ain'                  Aeq';
  //    -diag(lambda)*Ain -diag(fx)             zeros(m, size(Aeq,1));
  //    Aeq               zeros(size(Aeq,1), m) zeros(size(Aeq,1), size(Aeq,1))]
  // This implementation expects that Ain and Aeq are stored in the first
  // block-row. The first and third rows are automatically set when the user
  // sets Aeq and Ain. The second row is calculated below, performing more
  // operations, but saving memory. Only the second row will change.

  // Build Second Block-Row of F
  // Useful blocks
  const auto& x = y_->block(0, 0, n_, 1);  // Primal variables
  const auto& lambdas = y_->block(n_, 0, m_, 1);  // Dual variables (ineq)
  const auto& nu = y_->block(n_ + m_, 0, p_, 1);  // Dual variables (eq)
  const auto& Ain_t = F_->block(0, n_, n_, m_);  // Ain transpose
  const auto& Aeq_t = F_->block(0, n_ + m_, n_, p_);  // Aeq transpose

  Matrix<Scalar, Dynamic, 1> fx = Ain_t.transpose()*x - *bin_;
  DiagonalMatrix<Scalar, Dynamic, Dynamic> diag_lambdas(lambdas);

  // F_->block(n_, 0, m_, n_) = -1.0*diag_lambdas*Ain_t.transpose();
  // Copying block manually as Eigen does not support block assignment to Sparse
  // Matrices.
  Eigen::SparseMatrix<Scalar> first_block = -1.0*diag_lambdas*Ain_t.transpose();
  for (int i = 0; i < m_; i++) {
    for (int j = 0; j < n_; j++) {
      F_->coeffRef(n_ + i, j) = first_block.coeff(i, j);
    }
    F_->coeffRef(n_ + i, n_ + i) = -fx(i);
  }
  // Build Residuals
  *eta = static_cast<Scalar>(-fx.dot(lambdas));
  *t = mu * m_ / *eta;
  Scalar t_inv = static_cast<Scalar>(1.0 / *t);
  calculateResiduals(*y_, diag_lambdas, fx, t_inv, residuals);
  *rd_norm = residuals->block(0, 0, n_, 1).norm();
  *rp_norm = residuals->block(n_ + m_, 0, p_, 1).norm();
}

template <typename Scalar>
void
SparsePrimalDualLP<Scalar>::calculateResiduals(
    const Matrix<Scalar, Dynamic, 1>& y,
    const DiagonalMatrix<Scalar, Dynamic, Dynamic>& diag,
    const Matrix<Scalar, Dynamic, 1>& fx,
    const Scalar t_inv,
    Matrix<Scalar, Dynamic, 1>* residuals) {
  // Useful blocks
  const auto& x = y.block(0, 0, n_, 1);  // Primal variables
  const auto& lambdas = y.block(n_, 0, m_, 1);  // Dual variables (ineq)
  const auto& nu = y.block(n_ + m_, 0, p_, 1);  // Dual variables (eq)
  const auto& Ain_t = F_->block(0, n_, n_, m_);  // Ain transpose
  const auto& Aeq_t = F_->block(0, n_ + m_, n_, p_);  // Aeq transpose

  // Build Residuals
  // Dual residuals
  residuals->block(0, 0, n_, 1) = F_->block(0, 0, n_, l_)*y + *c_;
  // Central residuals
  residuals->block(n_, 0, m_, 1) =
       -t_inv*Matrix<Scalar, Dynamic, 1>::Ones(m_) - diag*fx;
  residuals->block(n_ + m_, 0, p_, 1) =
      Aeq_t.transpose()*y.block(0, 0, n_, 1) - *beq_;  // Primal residuals
}

template <typename Scalar>
bool
SparsePrimalDualLP<Scalar>::backTracking(
    const Scalar t_inv,
    const Scalar r_norm,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& ynt,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* y_plus,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* r_plus,
    Scalar* s) {
  // Find s_max
  // See page 613 from Boyd's book on Convex optimization
  *s = static_cast<Scalar>(1.00);
  const Scalar* ynt_data = ynt.data();
  for (uint i = 0; i < m_; i++) {
    if (ynt_data[n_ + i] < static_cast<Scalar>(0.0)) {
      Scalar val = -(*y_)(n_ + i) / ynt_data[n_ + i];
      if (val < *s) *s = val;
    }
  }

  // Do regular backtracking
  const Scalar zero = static_cast<Scalar>(0.0);
  *s *= static_cast<Scalar>(0.99);
  Scalar rhs;
  Scalar lhs;
  uint ninf = 0;  // Number of infeasible ineq.

  // Useful blocks
  const auto& x_vec = y_plus->block(0, 0, n_, 1);  // x plus
  const auto& lambdas_vec = y_plus->block(n_, 0, m_, 1);  // Lambdas plus
  const auto& Ain_t = F_->block(0, n_, n_, m_);  // Ain transpose

  do {
    ninf = 0;
    *y_plus = (*y_) + (*s) * ynt;
    // Calculate Residuals with y_plus
    Matrix<Scalar, Dynamic, 1> fx = Ain_t.transpose()*x_vec - *bin_;
    DiagonalMatrix<Scalar, Dynamic, Dynamic> diag_lambdas(lambdas_vec);
    calculateResiduals(*y_plus, diag_lambdas, fx, t_inv, r_plus);
    // Calculate lhs and rhs
    rhs = static_cast<Scalar>(1.0 - this->options.alpha_ * (*s)) * r_norm;
    lhs = r_plus->norm();
    // Calculate number of infeasible inequalities
    for (int k = 0; k < m_; k++) if (fx(k) >= 0) ninf++;
    // Calculate new s
    *s *= this->options.beta_;
  } while (lhs > rhs || ninf > 0);
  // Since we always do an extra multiplication (CHANGE)
  *s /= this->options.beta_;
  return true;
}
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_SPARSE_PRIMAL_DUAL_LP_IMPL_H_
