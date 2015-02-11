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

#ifndef OPTIMO_SOLVERS_NEWTON_IMPL_H_
#define OPTIMO_SOLVERS_NEWTON_IMPL_H_

namespace optimo {
namespace solvers {
//////////////////////////////////////////
// Newton Methods Implementation        //
template <typename Scalar, uint m, uint n>
TERMINATION_TYPE
Newton<Scalar, m, n>::operator()(const Problem<Scalar, m, n>& problem,
                                 Eigen::Matrix<Scalar, n, 1>* x,
                                 Scalar* min_value) {
  uint iter = 0;
  Eigen::Matrix<Scalar, n, n> hessian = Eigen::Matrix<Scalar, n, n>::Zero();
  Eigen::Matrix<Scalar, n, 1> gradient = Eigen::Matrix<Scalar, n, 1>::Zero();
  const GradientFunctor<Scalar, n>& gradient_functor = problem.gradient;
  const HessianFunctor<Scalar, n>& hessian_functor = problem.hessian;
  const ObjectiveFunctor<Scalar, n>& objective_functor = problem.objective;
  //  TERMINATION_TYPE r;

  while (++iter < this->options.max_iter_) {
    gradient_functor(*x, &gradient);
    hessian_functor(*x, &hessian);
    auto r = (*method_)(problem, hessian, gradient,
                        this->options.epsilon_, x, min_value);
    if (r == SOLVED) return r;
    else if (r == FAIL_NAN_INF) return r;
  }

  return MAX_ITERATIONS;
}

template <typename Scalar, uint m, uint n>
TERMINATION_TYPE
Newton<Scalar, m, n>::UnconstrainedMethod::operator()(
    const Problem<Scalar, m, n>& problem,
    const Eigen::Matrix<Scalar, n, n>& hessian,
    const Eigen::Matrix<Scalar, n, 1>& gradient,
    const Scalar& epsilon,
    Eigen::Matrix<Scalar, n, 1>* x,
    Scalar* min_value) {

  // 1. Compute Newton Step (Depending on the problem and type)
  const ObjectiveFunctor<Scalar, n>& objective_functor = problem.objective;
  this->xnt_ = -hessian.ldlt().solve(gradient);  // Newton Step
  lambda_sqrd_ = static_cast<Scalar>(-0.5*gradient.transpose()*this->xnt_);

  // 2. Stop if NAN or INF or 0.5*lambda_sqrd <= epsilon
  if (!optimo::utils::isFinite<Scalar, n, 1>(gradient) ||
      !optimo::utils::isFinite<Scalar, n, n>(hessian)) {
    return FAIL_NAN_INF;
  }
  // TODO(vfragoso): Check if it should be AND instead of an OR
  if (lambda_sqrd_ <= epsilon &&
      gradient.norm() < epsilon) {
    *min_value = objective_functor(*x);
    return SOLVED;
  }

  // 3. Line Search
  t_ = line_search(objective_functor, *x, gradient);

  // 4. Update x
  *x += t_*this->xnt_;
  return NOT_SOLVED;
}

template <typename Scalar, uint m, uint n>
Scalar Newton<Scalar, m, n>::UnconstrainedMethod::line_search(
    const ObjectiveFunctor<Scalar, n>& objective,
    const Eigen::Matrix<Scalar, n, 1>& x,
    const Eigen::Matrix<Scalar, n, 1>& g) {
  t_ = static_cast<Scalar>(1.0);
  fx_ = objective(x);
  Scalar g_xnt = static_cast<Scalar>(-0.5)*lambda_sqrd_;
  Scalar g_xnt_scalar = this->alpha_*g_xnt;
  Scalar lhs, rhs;
  do {
    lhs = objective(x + t_*this->xnt_);
    rhs = fx_ + t_*g_xnt_scalar;
    t_ *= this->beta_;
  } while (lhs > rhs);
  return t_;
}

template <typename Scalar, uint m, uint n>
TERMINATION_TYPE
Newton<Scalar, m, n>::EqualityConstrainedMethod::operator()(
    const Problem<Scalar, m, n>& problem,
    const Eigen::Matrix<Scalar, n, n>& hessian,
    const Eigen::Matrix<Scalar, n, 1>& gradient,
    const Scalar& epsilon,
    Eigen::Matrix<Scalar, n, 1>* x,
    Scalar* min_value) {
  const ObjectiveFunctor<Scalar, n>& objective_functor = problem.objective;
  const Eigen::Matrix<Scalar, m, n>& A = problem.A;  // Eq. Constraints
  const Eigen::Matrix<Scalar, m, 1>& b = problem.b;  // Eq. Constraints

  // 1. Compute Newton Step (Depending on the problem and type)
  const uint l = n + m;
  Eigen::Matrix<Scalar, l, l> AA;
  Eigen::Matrix<Scalar, l, 1> bb;
  AA.template block<n, n>(0, 0) = hessian;
  AA.template block<m, n>(n, 0) = A;
  AA.template block<n, m>(0, n) = A.transpose();
  AA.template block<m, m>(n, n) = Eigen::Matrix<Scalar, m, m>::Zero();
  bb.template block<n, 1>(0, 0) = -gradient;
  bb.template block<m, 1>(n, 0) = Eigen::Matrix<Scalar, m, 1>::Zero();
  Eigen::Matrix<Scalar, l, 1> nt = AA.ldlt().solve(bb);
  this->xnt_ = nt.template block<n, 1>(0, 0);

  // 2. Stop if NAN or INF or 0.5*lambda_sqrd <= epsilon or Gradient <= 0
  if (!optimo::utils::isFinite<Scalar, n, 1>(gradient) ||
      !optimo::utils::isFinite<Scalar, n, n>(hessian)) {
    *min_value = objective_functor(*x);
    return FAIL_NAN_INF;
  }
  this->lambda_sqrd_ = this->xnt_.transpose()*hessian*this->xnt_;
  if (0.5*this->lambda_sqrd_ <= epsilon) {
    *min_value = objective_functor(*x);
    return SOLVED;
  }

  // 3. Line Search
  this->t_ = this->line_search(objective_functor, *x, gradient);

  // 4. Update x
  *x += this->t_*this->xnt_;
  return NOT_SOLVED;
}

template <typename Scalar, uint m, uint n>
Scalar Newton<Scalar, m, n>::InfeasibleMethod::line_search(
    const ObjectiveFunctor<Scalar, n>& objective,
    const Eigen::Matrix<Scalar, n, 1>& x,
    const Eigen::Matrix<Scalar, m, 1>& v,
    const Eigen::Matrix<Scalar, m, n>& A,
    const Eigen::Matrix<Scalar, m, 1>& b,
    const Eigen::Matrix<Scalar, n, 1>& g) {
  const uint l = n + m;
  this->t_ = static_cast<Scalar>(1.0);
  Eigen::Matrix<Scalar, l, 1> r;
  Eigen::Matrix<Scalar, l, 1> rnt;
  r.template block<n, 1>(0, 0) = g + A.transpose()*v;
  r.template block<m, 1>(n, 0) = A*x - b;

  do {
    rnt.template block<n, 1>(0, 0) =
        r.template block<n, 1>(0, 0) + this->t_*A.transpose()*this->vnt_;
    rnt.template block<m, 1>(n, 0) =
        -r.template block<m, 1>(n, 0) - this->t_*A*this->xnt_;
    this->t_ *= this->beta_;
  } while (rnt.norm() > (1 - this->alpha_*t_)*r.norm());

  return t_;
}

template <typename Scalar, uint m, uint n>
TERMINATION_TYPE
Newton<Scalar, m, n>::InfeasibleMethod::operator()(
    const Problem<Scalar, m, n>& problem,
    const Eigen::Matrix<Scalar, n, n>& hessian,
    const Eigen::Matrix<Scalar, n, 1>& gradient,
    const Scalar& epsilon,
    Eigen::Matrix<Scalar, n, 1>* x,
    Scalar* min_value) {
  const ObjectiveFunctor<Scalar, n>& objective_functor = problem.objective;
  const Eigen::Matrix<Scalar, m, n>& A = problem.A;  // Eq. Constraints
  const Eigen::Matrix<Scalar, m, 1>& b = problem.b;  // Eq. Constraints

  // 1. Compute Newton Step (Depending on the problem and type)
  const uint l = n + m;
  Eigen::Matrix<Scalar, l, l> AA;
  Eigen::Matrix<Scalar, l, 1> bb;
  AA.template block<n, n>(0, 0) = hessian;
  AA.template block<m, n>(n, 0) = A;
  AA.template block<n, m>(0, n) = A.transpose();
  AA.template block<m, m>(n, n) = Eigen::Matrix<Scalar, m, m>::Zero();
  bb.template block<n, 1>(0, 0) = -gradient;
  bb.template block<m, 1>(n, 0) = b - A*(*x);
  Eigen::Matrix<Scalar, l, 1> nt = AA.ldlt().solve(bb);
  this->xnt_ = nt.template block<n, 1>(0, 0);
  this->vnt_ = nt.template block<m, 1>(n, 0) - (*v_);

  // 2. Stop if NAN or INF or 0.5*lambda_sqrd <= epsilon
  if (!optimo::utils::isFinite<Scalar, n, 1>(gradient) ||
      !optimo::utils::isFinite<Scalar, n, n>(hessian)) {
    return FAIL_NAN_INF;
  }

  // 3. Line Search
  this->t_ = line_search(objective_functor, *x, *v_, A, b, gradient);

  // 4. Update x
  *x += this->t_*this->xnt_;
  *v_ += this->t_*this->vnt_;

  // 5. Stop?
  Eigen::Matrix<Scalar, l, 1> r;
  r.template block<n, 1>(0, 0) = gradient + A.transpose()*(*v_);
  Eigen::Matrix<Scalar, m, 1> rprimal = A*(*x) - b;
  r.template block<m, 1>(n, 0) = rprimal;

  if (rprimal.norm() <= epsilon && r.norm() <= epsilon) {
    *min_value = objective_functor(*x);
    return SOLVED;
  }

  return NOT_SOLVED;
}
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_NEWTON_IMPL_H_
