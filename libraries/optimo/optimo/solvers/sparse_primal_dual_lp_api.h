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

#ifndef OPTIMO_SOLVERS_SPARSE_PRIMAL_DUAL_LP_API_H_
#define OPTIMO_SOLVERS_SPARSE_PRIMAL_DUAL_LP_API_H_

#include <glog/logging.h>
#include <algorithm>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SPQRSupport>
#include "optimo/solvers/solver.h"

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::DiagonalMatrix;
using Eigen::Dynamic;
using Eigen::SparseMatrix;

// Sparse Version
/// Solves a constrained linear program (LP)

/// Solves the following constrained Linear Program (LP):
/// \f{eqnarray*}{
///    \underset{\mathbf{x}}{\text{minimize}} & \mathbf{c}^T \mathbf{x} & \\
///    \text{subject to} & & \\
///    A_{\text{in}} \mathbf{x} & \leq & \mathbf{b}_{\text{in}} \\
///    A_{\text{eq}} \mathbf{x} & = & \mathbf{b}_{\text{eq}} \\
/// \f}
/// where \f$A_{\text{in}} \in \mathbb{R}^{m \times n} \f$,
/// \f$A_{\text{eq}} \in \mathbb{R}^{p \times n} \f$,
/// \f$\mathbf{b}_{\text{in}} \in \mathbb{R}^{m}\f$,
/// \f$\mathbf{b}_{\text{eq}} \in \mathbb{R}^{p}\f$,
/// \f$\mathbf{c} \in \mathbb{R}^n\f$, and \f$ \mathbf{x} \in \mathbb{R}^n \f$.
/// 
/// This LP solver implements a primal dual interior point method where
/// - n: Number of unkowns
/// - m: Number of inequalities
/// - p: Number of equalities
/// 
/// The solver requires a vector
/// \f[
///  \mathbf{y} =
///  \begin{bmatrix}
///  \mathbf{x}^T & \lambda^T & \nu^T
///  \end{bmatrix}^T
/// \f]
/// where the primal and dual variables are stored.
template <typename Scalar = double>
class SparsePrimalDualLP : public Solver<Scalar> {
 public:
  /// Sparse LP params
  class Params {
    friend SparsePrimalDualLP;
   public:
    /// LP parameters constructor
    Params(const uint n  /**< Number of unkowns */,
           const uint m /**< Number of inequality constraints */,
           const uint p /**< Number of equality constraints */) :
        n_(n), m_(m), p_(p), l_(n + m + p),
        beq(p), bin(m), c(n), F_(n + m + p, n + m + p) {
      bin.setConstant(0.0);
      beq.setConstant(0.0);
      c.setConstant(0.0);
    }

    // Destructor
    virtual ~Params(void) { }

    /// This function will reserve memory for the sparse matrix (optional).
    inline void reserve(void) {
      const int nnz = 2*(m_ + p_)*n_ + m_*m_;  // Estimated nnz elements
      F_.reserve(nnz);  // Reserve memory
    }

    /// Reserve memory per column for the sparse matrix (optional).
    inline bool reserve(const std::vector<uint>& col_sizes) {
      if (col_sizes.size() != l_) return false;
      F_.reserve(col_sizes);
      return true;
    }

    // Public members
    Matrix<Scalar, Dynamic, 1> beq;  ///< Vector for equalities
    Matrix<Scalar, Dynamic, 1> bin;  ///< Vector for inequalities
    Matrix<Scalar, Dynamic, 1> c;  ///< Cost vector of the LP

    // Placing matrices Ain, Aeq
    /// Method that inserts an element to matrix \f$A_{\text{eq}}\f$
    inline bool setAeqElement(const uint i, const uint j, const Scalar& val) {
      if (i >= p_ || j >= n_) return false;
      if (val == static_cast<Scalar>(0.0)) return false;
      const uint f_col = n_ + m_ + i;  // It is transposed in F
      const uint f_row = j;  // It is transposed
      F_.coeffRef(f_col, f_row) = val;  // To first block-row of F
      F_.coeffRef(f_row, f_col) = val;  // To third block-row of F
      return true;
    }

    /// Method that inserts an element to matrix \f$A_{\text{in}}\f$
    inline bool setAinElement(const uint i, const uint j, const Scalar& val) {
      if (i >= this->m_ || j >= this->n_) return false;
      if (val == static_cast<Scalar>(0.0)) return false;
      const uint f_col = n_ + i;
      const uint f_row = j;
      F_.coeffRef(f_row, f_col) = val;
      return true;
    }

   private:
    SparseMatrix<Scalar> F_;  // KKT system
    const uint n_;  // Num. of variables
    const uint m_;  // Num. of ineq. constraints
    const uint p_;  // Nun. of eq. constraints
    const uint l_;  // l = n + m + p
  };

  // Constructor
  SparsePrimalDualLP(void) { }

  // Destructor
  virtual ~SparsePrimalDualLP(void) {
    F_ = nullptr;
    beq_ = nullptr;
    bin_ = nullptr;
    c_ = nullptr;
    y_ = nullptr;
  }

  /// Solves the LP problem
  TERMINATION_TYPE
  operator()(Params* params  /**< LP parameters*/,
             Matrix<Scalar, Dynamic, 1>* y /**< Vector of unkowns*/,
             Scalar* min_value /**< Minimum value*/);

 protected:
  inline bool feasibleStartingPoint(const Matrix<Scalar, Dynamic, 1>& y) {
    // Useful blocks
    const auto& x = y.block(0, 0, n_, 1);  // Primal variables
    const auto& Ain_t = F_->block(0, n_, n_, m_);  // Ain transpose
    const auto& Aeq_t = F_->block(0, n_ + m_, n_, p_);  // Aeq transpose
    Matrix<Scalar, Dynamic, 1> ineq_x = Ain_t.transpose()*x - *bin_;
    for (int i = 0; i < m_; i++) if (ineq_x(i) > 0.0) return false;
    // // TODO(vfragoso): Check if we can ignore the eq. constraints
    // Matrix<Scalar, Dynamic, 1> eq_x = Aeq_t.transpose()*x - *beq_;
    // for (int i = 0; i < p_; i++) {
    //   if (fabs(eq_x(i)) > this->options.eps_feas_) {
    //     return false;
    //   }
    // }
    return true;
  }

  bool backTracking(const Scalar t,
                    const Scalar r_norm,
                    const Matrix<Scalar, Dynamic, 1>& ynt,
                    Matrix<Scalar, Dynamic, 1>* y_plus,
                    Matrix<Scalar, Dynamic, 1>* r_plus,
                    Scalar* s);

  void calculateResiduals(
      const Matrix<Scalar, Dynamic, 1>& y,
      const DiagonalMatrix<Scalar, Dynamic, Dynamic>& diag,
      const Matrix<Scalar, Dynamic, 1> & fx,
      const Scalar t_inv,
      Matrix<Scalar, Dynamic, 1>* residuals);

  void buildKKTSystem(const double mu,
                      double* eta,
                      double* t,
                      Matrix<Scalar, Dynamic, 1>* residuals,
                      double* rd_norm,
                      double* rp_norm);

  inline TERMINATION_TYPE
  mapResult(const Eigen::SPQR<Eigen::SparseMatrix<Scalar> >& solver) {
    switch (solver.info()) {
      case Eigen::NumericalIssue:
        return NUMERICAL_ISSUE;
      case Eigen::NoConvergence:
        return NO_CONVERGENCE;
      case Eigen::InvalidInput:
        return INVALID_ARGUMENTS;
      case Eigen::Success:
        return SOLVED;
    }
  }

  inline bool isF_empty(void) {
    Scalar acc = static_cast<Scalar>(0.0);
    for (int i = 0; i < F_->outerSize(); i++) {
      for (typename SparseMatrix<Scalar>::InnerIterator it(*F_, i); it; ++it) {
        acc += static_cast<Scalar>(fabs(it.value()));
      }
    }
    return acc == static_cast<Scalar>(0.0);
  }

 private:
  uint n_;  // Num. of variables
  uint m_;  // Num. of ineq. constraints
  uint p_;  // Nun. of eq. constraints
  uint l_;  // l = n + m + p

  // Members so that the user can fill out the matrices
  // J(x) = c'*x
  SparseMatrix<Scalar>* F_ = nullptr;  // KKT system
  Matrix<Scalar, Dynamic, 1>* beq_ = nullptr;  // Vector of equalities
  Matrix<Scalar, Dynamic, 1>* bin_ = nullptr;  // Vector of equalities
  Matrix<Scalar, Dynamic, 1>* c_ = nullptr;  // Cost vector
  Matrix<Scalar, Dynamic, 1>* y_ = nullptr;  // Unknowns [x lambdas nu]'
};
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_SPARSE_PRIMAL_DUAL_LP_API_H_
