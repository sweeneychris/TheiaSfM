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

#ifndef OPTIMO_SOLVERS_BFGS_IMPL_H_
#define OPTIMO_SOLVERS_BFGS_IMPL_H_

#include "optimo/utils/matrix_utils.h"

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::Dynamic;

// BFGS Solver Implementation using Eigen Dynamic matrices
template <typename Scalar>
TERMINATION_TYPE
BFGS<Scalar>::operator()(const ProblemLS<Scalar>& problem,
                         Matrix<Scalar, Dynamic, 1>* x,
                         Scalar* min_value) {
  if (!x || !min_value) return NOT_SOLVED;
  int iter = 0;
  const int n = x->rows();
  Matrix<Scalar, Dynamic, Dynamic> H(n, n);
  Matrix<Scalar, Dynamic, 1> gradient(n), past_gradient(n), y(n), s(n), p(n);
  const GradientFunctorLS<Scalar>& gradient_functor = problem.gradient;
  const ObjectiveFunctorLS<Scalar>& objective_functor = problem.objective;
  TERMINATION_TYPE r;
  Scalar t = static_cast<Scalar>(1.0);  // Scalar for line search

  // Initialize gradient and ''Hessian''
  gradient.setZero();
  gradient_functor(*x, &gradient);  // Initialize gradient
  H.setIdentity();  // Initialize inverse of ''Hessian''

  while (++iter < this->options.max_iter_) {
    // Check if we converged
    if (!optimo::utils::isFinite<Scalar>(gradient) ||
        !optimo::utils::isFinite<Scalar>(H)) {
      return FAIL_NAN_INF;
    }
    // If gradient's norm is small or cost did not change much
    if (gradient.norm() < this->options.epsilon_) {
      *min_value = objective_functor(*x);
      return SOLVED;
    }

    // 1. Compute descent direction
    p = -H*gradient;

    // 1. Compute line search
    t = line_search(objective_functor, *x, p, gradient);

    // 2. Update x
    s = t*p;
    *x += s;

    // Update gradient and ''Hessian''
    past_gradient = gradient;  // Copy last gradient
    gradient_functor(*x, &gradient);  // Update gradient
    y = gradient - past_gradient;

    // Update Hessian via a secant method
    // Note: We use dot() instead of s.transpose()*y + y.transpose()*H*y as
    // Eigen struggled to evaluate that expression correctly with GNU GCC.
    Scalar den = y.transpose()*s;
    Scalar den_sqrd = den*den;
    Scalar scalar1 = (s.dot(y) + y.dot(H*y)) / den_sqrd;
    H += scalar1*s*s.transpose() -
        (H*y*s.transpose() + s*y.transpose()*H) / den;
  }

  return MAX_ITERATIONS;
}

// Line search: Implementing Backtracking as described in
// Boyd and Vandenberghe's book on Convex Optimization (pg 464)
template <typename Scalar>
Scalar BFGS<Scalar>::line_search(
    const ObjectiveFunctorLS<Scalar>& objective,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& p,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& g) {
  Scalar t = static_cast<Scalar>(1.0);
  Scalar fx = objective(x);
  Scalar lhs, rhs;
  Scalar update_scalar = this->options.alpha_*g.transpose()*p;
  do {
    lhs = objective(x + t*p);
    rhs = fx + t*update_scalar;
    t *= this->options.beta_;
  } while (lhs > rhs);
  return t;
}
}  // solvers
}  // optimo

#endif  // OPTIMO_SOLVERS_BFGS_IMPL_H_
