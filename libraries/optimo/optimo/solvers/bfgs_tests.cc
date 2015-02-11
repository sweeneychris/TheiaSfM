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


#include "optimo/solvers/bfgs.h"
#include <glog/logging.h>
#include "gtest/gtest.h"

namespace optimo {
namespace solvers {
namespace {
using Eigen::Matrix;
using Eigen::Dynamic;
// Quadratic Objective
template <typename Scalar, uint n>
class QuadraticObjective : public ObjectiveFunctorLS<Scalar> {
 public:
  QuadraticObjective(const Matrix<Scalar, n, n>& _Q,
                     const Matrix<Scalar, n, 1>& _P,
                     const Scalar& _c) {
    Q = _Q;
    P = _P;
    c = _c;
  }

  ~QuadraticObjective(void) { }

  virtual Scalar operator()(const Matrix<Scalar, Dynamic, 1>& x) const {
    Matrix<Scalar, 1, 1> sum = x.transpose()*Q*x + P.transpose()*x;
    return sum(0, 0) + c;
  }

 protected:
  Matrix<Scalar, n, n> Q;
  Matrix<Scalar, n, 1> P;
  Scalar c;
};

// Quadratic Program
template <typename Scalar, uint n>
class QuadraticProgram : public ProblemLS<Scalar> {
 public:
  QuadraticProgram(const QuadraticObjective<Scalar, n>& qobj,
                   const GradientFunctorLS<Scalar>& g,
                   const HessianFunctorLS<Scalar>& h) :
      ProblemLS<Scalar>(qobj, g, h) {
  }

  virtual ~QuadraticProgram(void) { }
};

// Gradient Functor
template <typename Scalar, uint n>
class QGradient : public GradientFunctorLS<Scalar> {
 public:
  QGradient(const Matrix<Scalar, n, n>& _Q,
            const Matrix<Scalar, n, 1>& _P) {
    Q = _Q;
    P = _P;
  }

  virtual void operator()(const Matrix<Scalar, Dynamic, 1>& x,
                          Matrix<Scalar, Dynamic, 1>* g) const {
    *g = 2*Q*x + P;
  }

 protected:
  Matrix<Scalar, n, n> Q;
  Matrix<Scalar, n, 1> P;
};

// Hessian Functor
template <typename Scalar, uint n>
class QHessian : public HessianFunctorLS<Scalar> {
 public:
  QHessian(const Matrix<Scalar, n, n>& _Q) {
    Q = _Q;
  }

  virtual void operator()(const Matrix<Scalar, Dynamic, 1>& x,
                          Matrix<Scalar, Dynamic, Dynamic>* h) const {
    *h = 2*Q;
  }

 protected:
  Matrix<Scalar, n, n> Q;
};

// Rosenbrock problem
template <typename Scalar>
class RosenbrockObjective : public ObjectiveFunctorLS<Scalar> {
 public:
  virtual Scalar operator()(const Matrix<Scalar, Dynamic, 1>& x) const {
    const Scalar term1 = 1 - x(0);
    const Scalar term2 = x(1) - x(0)*x(0);
    return term1*term1 + static_cast<Scalar>(100)*term2*term2;
  }
};

template <typename Scalar>
class RosenbrockGradient : public GradientFunctorLS<Scalar> {
 public:
  virtual void operator()(const Matrix<Scalar, Dynamic, 1>& x,
                          Matrix<Scalar, Dynamic, 1>* g) const {
    (*g)(0) = -2*(1 - x(0)) - 400*(x(1) - x(0)*x(0))*x(0);
    (*g)(1) = 200*(x(1) - x(0)*x(0));
  }
};

template <typename Scalar>
class RosenbrockHessian : public HessianFunctorLS<Scalar> {
 public:
  virtual void operator()(const Matrix<Scalar, Dynamic, 1>& x,
                          Matrix<Scalar, Dynamic, Dynamic>* h) const {
  }
};

template <typename Scalar, uint n = 2>
class RosenbrockProblem : public ProblemLS<Scalar> {
 public:
  RosenbrockProblem(const RosenbrockObjective<Scalar>& obj,
                    const RosenbrockGradient<Scalar>& g,
                    const RosenbrockHessian<Scalar>& h) :
      ProblemLS<Scalar>(obj, g, h) { }

  virtual ~RosenbrockProblem(void) { }
};
}  // namespace

TEST(BFGS, Quadratic_Program) {
  const int n = 3;  // Number of unkowns
  const uint m = 1;  // Number of constraints
  Matrix<double, n, n> Q;
  Matrix<double, n, 1> P;
  double c = 15.0;
  Matrix<double, Dynamic, 1> x(n);

  Q << 5, 0, 0,
       0, 4, 0,
       0, 0, 1.5;
  P << 0.01, 0.5, -0.1;
  x << 100, 100, 100;

  QuadraticObjective<double, n> qobj(Q, P, c);
  QGradient<double, n> qgradient(Q, P);
  QHessian<double, n> qhessian(Q);

  QuadraticProgram<double, n> qprogram(qobj, qgradient, qhessian);

  const double epsilon = 1e-6;
  const int max_iter = 50;
  const double alpha = 0.01;
  const double beta = 0.5;

  BFGS<double> bfgs;
  bfgs.options.epsilon_ = epsilon;
  bfgs.options.max_iter_ = max_iter;
  bfgs.options.alpha_ = alpha;
  bfgs.options.beta_ =beta;
  double min_val;
  auto res = bfgs(qprogram, &x, &min_val);
  VLOG(1) << "Min_val: " << min_val
          << " Minimizer: [" << x.transpose() << "]"
          << " Result: " << res;
  ASSERT_EQ(res, 0) << "Algorithm did not converge: res= " << res;
  ASSERT_LT(min_val - c, 0.1);
}

TEST(BFGS, Rosenbrock) {
  RosenbrockObjective<double> obj;
  RosenbrockGradient<double> g;
  RosenbrockHessian<double> h;
  RosenbrockProblem<double> problem(obj, g, h);

  const double epsilon = 1e-6;
  const int max_iter = 100;
  const double alpha = 0.01;
  const double beta = 0.5;
  BFGS<double> bfgs;
  bfgs.options.epsilon_ = epsilon;
  bfgs.options.max_iter_ = max_iter;
  bfgs.options.alpha_ = alpha;
  bfgs.options.beta_ =beta;
  double min_val;
  Matrix<double, Dynamic, 1> x(2);
  x(0) = 0.5;
  x(1) = 3.0;
  auto res = bfgs(problem, &x, &min_val);
  VLOG(1) << "Min_val: " << min_val
          << " Minimizer: [" << x.transpose() << "]"
          << " Result: " << res;
  ASSERT_EQ(res, 0) << "Algorithm did not converge: res= " << res;
  ASSERT_NEAR(0.0, min_val, 1e-3);
  ASSERT_NEAR(x(0), 1.0, 1e-3);
  ASSERT_NEAR(x(1), 1.0, 1e-3);
}
}  // solvers
}  // optimo
