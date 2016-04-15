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

#ifndef OPTIMO_SOLVERS_PRIMAL_DUAL_QP_API_H_
#define OPTIMO_SOLVERS_PRIMAL_DUAL_QP_API_H_

#include <Eigen/Core>
#include <glog/logging.h>
#include "optimo/solvers/solver.h"

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::Dynamic;
// Solves a constrained quadratic program (QP).

// Solves the following constrained quadratic program (QP):
// \f{eqnarray*}{
// \underset{\mathbf{x}}{\text{minimize}} &
//  \frac{\mathbf{x}^T Q \mathbf{x}}{2} + & \mathbf{d}^T\mathbf{x}
// \text{subject to} & &
//  A_{\text{eq}} \mathbf{x} & = & \mathbf{b}_{\text{eq}}
//  A_{\text{in}} \mathbf{x} & \leq & \mathbf{b}_{\text in}
// \f}
// where \f$A_{\text{in}} \in \mathbb{R}^{m \times n} \f$,
// \f$A_{\text{eq}} \in \mathbb{R}^{p \times n} \f$,
// \f$\mathbf{b}_{\text{in}} \in \mathbb{R}^{m}\f$,
// \f$\mathbf{b}_{\text{eq}} \in \mathbb{R}^{p}\f$,
// \f$\mathbf{d} \in \mathbb{R}^n\f$, \f$ Q \in \mathbb{R}^{n \times n} \f$
// and \f$ \mathbf{x} \in \mathbb{R}^n \f$.
//
// This QP solver implements a primal dual interior point method where
// - n: Number of unkowns
// - m: Number of inequalities
// - p: Number of equalities
// This solver is recomended for small problems (n + m + p < 100). For larger
// problems consider using the Sparse implementation (SparsePrimalDualQP).
//
// The solver requires a vector
// \f[
//  \mathbf{y} =
//  \begin{bmatrix}
//  \mathbf{x}^T & \lambda^T & \nu^T
//  \end{bmatrix}^T
// \f]
// where the primal and dual variables are stored.
template <typename Scalar, uint n, uint m, uint p>
class PrimalDualQP : public Solver<Scalar> {
 public:
  // QP parameters defining the problem
  struct Params {
    Params(void) : Q(n, n), d(n), Aeq(p, n), beq(p), Ain(m, n), bin(m) { }

    // QP objective params
    Matrix<Scalar, Dynamic, Dynamic> Q;  // Q matrix
    Matrix<Scalar, Dynamic, 1> d;  // d vector
    // Constraints
    Matrix<Scalar, Dynamic, Dynamic> Aeq;  // Equality constraints
    Matrix<Scalar, Dynamic, 1> beq;  // Equality constraints vector
    Matrix<Scalar, Dynamic, Dynamic> Ain;  // Inequality constraints
    Matrix<Scalar, Dynamic, 1> bin;  // Inequality constraints vector
  };

  // Constructor
  PrimalDualQP(void) : ones_(Matrix<Scalar, Dynamic, 1>::Ones(m)) { }

  // Destructor
  virtual ~PrimalDualQP(void) { }

  // Solve the QP program
  TERMINATION_TYPE
  operator()(const typename PrimalDualQP<Scalar, n, m, p>::Params& params,
             Matrix<Scalar, n, 1>* x,
             Scalar* min_value);

 protected:
  // Line search procedure (Backtracking)
  Scalar backTracking(const Scalar t_inv,
                      const Scalar r_norm,
                      const Matrix<Scalar, Dynamic, 1>& ynt,
                      const Matrix<Scalar, Dynamic, 1>& y,
                      const Matrix<Scalar, Dynamic, 1>& d,
                      const Matrix<Scalar, p, 1>& beq,
                      const Matrix<Scalar, m, 1>& bin,
                      Matrix<Scalar, Dynamic, 1>* y_plus,
                      Matrix<Scalar, Dynamic, 1>* r_plus,
                      bool* exit_flag);

  // Calculates residuals (primal and dual)
  inline
  void calculateResiduals(const Matrix<Scalar, Dynamic, 1>& y,
                          const Matrix<Scalar, Dynamic, 1>& d,
                          const Matrix<Scalar, Dynamic, 1>& diag,
                          const Matrix<Scalar, p, 1>& beq,
                          const Scalar t_inv,
                          Matrix<Scalar, Dynamic, 1>* r);

 private:
  Matrix<Scalar, Dynamic, Dynamic> F_;  // Matrix for KKT system
  Matrix<Scalar, Dynamic, 1> r_;  // Residual vector
  const Matrix<Scalar, Dynamic, 1> ones_;  // Helper vector
};
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_PRIMAL_DUAL_QP_API_H_
