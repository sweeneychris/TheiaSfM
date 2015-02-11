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

#ifndef OPTIMO_SOLVERS_GRADIENT_DESCENT_H_
#define OPTIMO_SOLVERS_GRADIENT_DESCENT_H_

#include <glog/logging.h>
#include <Eigen/Core>
#include "optimo/core/objects.h"
#include "optimo/core/objects_ls.h"
#include "optimo/utils/matrix_utils.h"
#include "optimo/solvers/solver.h"

namespace optimo {
namespace solvers {
/// Standard Gradient Descent method.

/// This class implements the following procedure:
/// 1. Computes a descent direction \f$\mathbf{g}_k=-\nabla f_0(\mathbf{x}_k)\f$
/// 2. Stop if the norm of the gradient is \f$\|\mathbf{g}_k\|\leq\epsilon\f$
/// 3. Take a step in the descent direction
/// \f[
/// \mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha\mathbf{g}_k
/// \f]
///
/// where \f$\alpha\f$ is the step size.
///
/// This implementation works well for small problems that fit in the
/// stack m < 16, n < 16. If the problem does not fit in the stack, consider
/// using the GradientDescentLS class.
template <typename Scalar, uint m, uint n>
class GradientDescent : public Solver<Scalar> {
 public:
  /// Constructor
  explicit GradientDescent(const Scalar step_size = 1e-3 /**< Step size*/) :
      step_size_(step_size) { }

  // Destructor
  virtual ~GradientDescent(void) { }

  /// Solve the problem
  TERMINATION_TYPE
  operator()(const Problem<Scalar, m, n>& problem,
             Eigen::Matrix<Scalar, n, 1>* x,
             Scalar* min_value);

 protected:
  const Scalar step_size_;  // Step size for update
};

template <typename Scalar, uint m, uint n>
TERMINATION_TYPE
GradientDescent<Scalar, m, n>::operator()(const Problem<Scalar, m, n>& problem,
                                          Eigen::Matrix<Scalar, n, 1>* x,
                                          Scalar* min_value) {
  uint iter = 0;
  Eigen::Matrix<Scalar, n, 1> gradient = Eigen::Matrix<Scalar, n, 1>::Zero();
  const GradientFunctor<Scalar, n>& gradient_functor = problem.gradient;
  const ObjectiveFunctor<Scalar, n>& objective_functor = problem.objective;
  TERMINATION_TYPE r;

  while (++iter < this->options.max_iter_) {
    // 1. Compute descent direction
    gradient_functor(*x, &gradient);

    // 1.1 Check that gradient does not have inf or nan elements.
    if (!optimo::utils::isFinite<Scalar, n, 1>(gradient)) {
      return FAIL_NAN_INF;
    }

    // 2. Check norm, if norm is really small then we are in a local optima
    if (gradient.norm() <= this->options.epsilon_) {
      // TODO(vfragoso): Check Hessian to see if it is a minima or maxima
      // if it is minima we are done, if not then the problem was not solved!
      r = SOLVED;
      *min_value = objective_functor(*x);
      return r;
    }

    // 2. Update
    // TODO(vfragoso): Implement a line-search method
    *x -= step_size_*gradient;
  }

  return MAX_ITERATIONS;
}

/// Standard Gradient Descent method.

/// This class implements the following procedure:
/// 1. Computes a descent direction \f$\mathbf{g}_k=-\nabla f_0(\mathbf{x}_k)\f$
/// 2. Stop if the norm of the gradient is \f$\|\mathbf{g}_k\|\leq\epsilon\f$
/// 3. Take a step in the descent direction
/// \f[
/// \mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha\mathbf{g}_k
/// \f]
///
/// where \f$\alpha\f$ is the step size.
template <typename Scalar>
class GradientDescentLS : public Solver<Scalar> {
 public:
  /// Constructor
  GradientDescentLS(const int n,  ///< Number of variables
                    const Scalar alpha = 1e-3  /**< Step size*/)
      : n_(n), step_size_(alpha) { }

  // Destructor
  virtual ~GradientDescentLS(void) { }

  /// Solve the problem
  TERMINATION_TYPE
  operator()(const ProblemLS<Scalar>& problem,
             Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* x,
             Scalar* min_value);

 protected:
  const int n_;  // Number of variables
  const Scalar step_size_;  // Step size for update
};

template <typename Scalar>
TERMINATION_TYPE
GradientDescentLS<Scalar>::operator()(
    const ProblemLS<Scalar>& problem,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* x,
    Scalar* min_value) {
  uint iter = 0;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> gradient(n_);
  gradient.setZero();
  const GradientFunctorLS<Scalar>& gradient_functor = problem.gradient;
  const ObjectiveFunctorLS<Scalar>& objective_functor = problem.objective;
  TERMINATION_TYPE r;

  while (++iter < this->options.max_iter_) {
    // 1. Compute descent direction
    gradient_functor(*x, &gradient);

    // 1.1 Check that gradient does not have inf or nan elements.
    if (!optimo::utils::isFinite<Scalar>(gradient)) {
      return FAIL_NAN_INF;
    }

    // 2. Check norm, if norm is really small then we are in a local optima
    if (gradient.norm() <= this->options.epsilon_) {
      r = SOLVED;
      *min_value = objective_functor(*x);
      return r;
    }

    // 2. Update
    // TODO(vfragoso): Implement a line-search method
    *x -= step_size_*gradient;
  }

  return MAX_ITERATIONS;
}
}  // solvers
}  // optimo

#endif  // OPTIMO_SOLVERS_GRADIENT_DESCENT_H_
