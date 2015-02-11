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

#ifndef OPTIMO_SOLVERS_BFGS_API_H_
#define OPTIMO_SOLVERS_BFGS_API_H_

#include "optimo/core/objects_ls.h"
#include <Eigen/Core>
#include "optimo/solvers/solver.h"

namespace optimo {
namespace solvers {
/// Implements the Broyden–Fletcher–Goldfarb–Shanno algorithm

/// Implementation of the Broyden–Fletcher–Goldfarb–Shanno quasi newton
/// algorithm. This algorithm solves unconstrained problems. The implementation
/// can be summarized as follows:
/// 1. Compute a descent direction by computing 
/// \f[
///  \mathbf{p}_k = -H_k \nabla f_0(\mathbf{x}_k)
/// \f]
/// 2. Compute a step size \f$ t_k \f$by a line search (Backtracking)
/// 3. Update: \f$ \mathbf{x}_{k+1} = \mathbf{x}_k + t_k \mathbf{p}_k \f$
/// 4. Set \f$ \mathbf{s}_k = t_k \mathbf{p}_k \f$
/// 5. Compute
/// \f[
/// \mathbf{y}_k = \nabla f_0(\mathbf{x}_{k+1}) - \nabla f_0(\mathbf{x}_k)
/// \f]
/// 4. Compute \f$ H_{k+1} \f$ as follows:
/// \f[
/// H_{k+1} = H_k + T_1 - T_2
/// \f]
/// where
/// \f{eqnarray*}{
/// T_1 &=& \alpha \frac{\mathbf{s}_k \mathbf{s}_k^T}{\beta^2} \\
/// T_2 &=&
///  \frac{H_k \mathbf{y}_k \mathbf{s}_k^T + \mathbf{s}_k \mathbf{y}_k^T H_k}
///   {\beta} \\
/// \alpha & = & \mathbf{s}_k^T \mathbf{y}_k + \mathbf{y}_k^T H_k\mathbf{y}_k \\
/// \beta & = & \mathbf{s}_k^T \mathbf{y} \\
/// \f}
/// and \f$ f_0 \f$ is the objective function. Further information can be found
/// in the Wikipedia 
/// <a href="http://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm">
/// article </a> of the BFGS this algoritm.
template <typename Scalar>
class BFGS : public Solver<Scalar> {
 public:
  // Constructor for BFGS Solver
  BFGS(void) { }

  // Destructor
  virtual ~BFGS(void) { }

  /// Solve the problem
  TERMINATION_TYPE
  operator()(const ProblemLS<Scalar>& problem,
             Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* x,
             Scalar* min_value);

 protected:
  /// Line search: Implements a Backtracking algorithm
  Scalar line_search(const ObjectiveFunctorLS<Scalar>& objective,
                     const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
                     const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& p,
                     const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& g);
};
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_BFGS_API_H_

