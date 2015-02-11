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

#ifndef OPTIMO_SOLVERS_PRIMAL_DUAL_LP_API_H_
#define OPTIMO_SOLVERS_PRIMAL_DUAL_LP_API_H_

#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include "optimo/solvers/solver.h"

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::DiagonalMatrix;
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
/// This solver is recomended for small problems (n + m + p < 100). For larger
/// problems consider using the SparsePrimalDualLP solver.
///
/// The solver requires a vector
/// \f[
///  \mathbf{y} =
///  \begin{bmatrix}
///  \mathbf{x}^T & \lambda^T & \nu^T
///  \end{bmatrix}^T
/// \f]
/// where the primal and dual variables are stored.
template <typename Scalar, uint n, uint m, uint p>
class PrimalDualLP : public Solver<Scalar> {
 public:
  /// Linear Program (LP) Parameters defining the problem
  struct Params {
    Params(void) : c(n), Ain(m, n), bin(m), Aeq(p, n), beq(p) {
      Ain.setConstant(static_cast<Scalar>(0.0));
      Aeq.setConstant(static_cast<Scalar>(0.0));
      c.setConstant(static_cast<Scalar>(0.0));
      bin.setConstant(static_cast<Scalar>(0.0));
      beq.setConstant(static_cast<Scalar>(0.0));
    }

    // LP Params
    Matrix<Scalar, Dynamic, 1> c;  ///< Costs vector
    Matrix<Scalar, Dynamic, Dynamic> Aeq;  ///< Equality constraints
    Matrix<Scalar, Dynamic, 1> beq;  ///< Equality constraints vector
    Matrix<Scalar, Dynamic, Dynamic> Ain;  ///< Inequality constraints
    Matrix<Scalar, Dynamic, 1> bin;  ///< Inequality constraints vector
  };

  // Constructor
  PrimalDualLP(void) { }

  // Destructor
  virtual ~PrimalDualLP(void) { }

  /// Solves the LP program defined in params.
  TERMINATION_TYPE
  operator()(const Params& params,
             Matrix<Scalar, Dynamic, 1>* y,
             Scalar* min_value);

 protected:
  // Calculate residuals
  void calculateResiduals(const Params& params,
                          const Matrix<Scalar, Dynamic, 1>& y,
                          const DiagonalMatrix<Scalar, Dynamic, Dynamic>& diag,
                          const Matrix<Scalar, Dynamic, 1>& fx,
                          const Scalar t_inv,
                          Matrix<Scalar, Dynamic, 1>* residuals);
  // Line Search
  Scalar backTracking(const double t_inv,
                      const double r_norm,
                      const Matrix<Scalar, Dynamic, 1>& ynt,
                      const Matrix<Scalar, Dynamic, 1>& y,
                      const Params& params);

  // Buld the KKT System
  void buildKKTSystem(const Params& params,
                      const Matrix<Scalar, Dynamic, 1>& y,
                      const Scalar mu,
                      Scalar* eta,
                      Scalar* t,
                      Matrix<Scalar, Dynamic, 1>* residuals,
                      Scalar* rd_norm,
                      Scalar* rp_norm);

 private:
  // Members
  Matrix<Scalar, Dynamic, Dynamic> F_;  // KKT system
  Matrix<Scalar, Dynamic, 1> r_;  // Residuals
  Matrix<Scalar, Dynamic, 1> y_plus_;  // Variables for line search
  Matrix<Scalar, Dynamic, 1> r_plus_;  // Residuals for line search
};
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_PRIMAL_DUAL_LP_API_H_
