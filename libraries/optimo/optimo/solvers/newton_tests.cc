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

#include "optimo/core/objects.h"
#include "optimo/solvers/newton.h"
#include "optimo/core/numerical_gradient.h"
#include "optimo/core/numerical_hessian.h"
#include <glog/logging.h>
#include "gtest/gtest.h"
#include <Eigen/LU>
#include <random>
#include <cmath>

namespace optimo {
namespace solvers {
namespace {
using Eigen::Matrix;
// Quadratic Objective
template <typename Scalar, uint n>
class QuadraticObjective : public ObjectiveFunctor<Scalar, n> {
 public:
  QuadraticObjective(const Matrix<Scalar, n, n>& _Q,
                     const Matrix<Scalar, n, 1>& _P,
                     const Scalar& _c) {
    Q = _Q;
    P = _P;
    c = _c;
  }

  ~QuadraticObjective(void) { }
  
  virtual Scalar operator()(const Matrix<Scalar, n, 1>& x) const {
    Matrix<Scalar, 1, 1> sum = x.transpose()*Q*x + P.transpose()*x;
    return sum(0,0) + c;
  }
  
 protected:
  Matrix<Scalar, n, n> Q;
  Matrix<Scalar, n, 1> P;
  Scalar c;
};

// Quadratic Program
template <typename Scalar, uint m, uint n>
class QuadraticProgram : public Problem<Scalar, m, n> {
 public:
  QuadraticProgram(const QuadraticObjective<Scalar, n>& qobj,
                   const NoConstraints<Scalar, m, n>& nc,
                   const Matrix<Scalar, m, n>& Aeq,
                   const Matrix<Scalar, m, 1>& beq,
                   const GradientFunctor<Scalar, n>& g,
                   const HessianFunctor<Scalar, n>& h) :
      Problem<Scalar, m, n>(qobj, nc, Aeq, beq, g, h) {
  }

  virtual ~QuadraticProgram(void) { }
};

// Gradient Functor
template <typename Scalar, uint n>
class QGradient : public GradientFunctor<Scalar, n> {
 public:
  QGradient(const Matrix<Scalar, n, n>& _Q,
            const Matrix<Scalar, n, 1>& _P) {
    Q = _Q;
    P = _P;
  }

  virtual void operator()(const Matrix<Scalar, n, 1>& x,
                          Matrix<Scalar, n, 1>* g) const {
    *g = 2*Q*x + P;
  }

 protected:
  Matrix<Scalar, n, n> Q;
  Matrix<Scalar, n, 1> P;
};

// Hessian Functor
template <typename Scalar, uint n>
class QHessian : public HessianFunctor<Scalar, n> {
 public:
  QHessian(const Matrix<Scalar, n, n>& _Q) {
    Q = _Q;
  }

  virtual void operator()(const Matrix<Scalar, n, 1>& x,
                          Matrix<Scalar, n, n>* h) const {
    *h = 2*Q;
  }

 protected:
  Matrix<Scalar, n, n> Q;
};

// Entropy Maximization Problem
template <typename Scalar, uint m, uint n>
struct EntropyObjective : public ObjectiveFunctor<Scalar, n> {
  virtual Scalar operator()(const Matrix<Scalar, n, 1>& x) const {
    Scalar f = static_cast<Scalar>(0.0);
    const Scalar* xdata = x.data();
    for (unsigned int i = 0; i < n; i++) {
      // TODO(vfragoso): Check if log returns NaN when 0 or x  < 0
      Scalar logx = static_cast<Scalar>(log(xdata[i]));
      if (std::isnan(logx)) return logx;
      f += xdata[i]*logx;
    }
    return f;
  }
};

// Gradient Functor
template <typename Scalar, uint n>
class EntropyGradient : public GradientFunctor<Scalar, n> {
  virtual void operator()(const Matrix<Scalar, n, 1>& x,
                          Matrix<Scalar, n, 1>* g) const {
    Scalar* gdata = g->data();
    const Scalar* xdata = x.data();
    for (unsigned int i = 0; i < n; i++) {
      gdata[i] = log(xdata[i]) + 1;
    }
  }
};

// Hessian Functor
template <typename Scalar, uint n>
class EntropyHessian : public HessianFunctor<Scalar, n> {
 public:
  virtual void operator()(const Matrix<Scalar, n, 1>& x,
                          Matrix<Scalar, n, n>* h) const {
    const Scalar* xdata = x.data();
    for (unsigned int i = 0; i < n; i++) {
      (*h)(i, i) = static_cast<Scalar>(1.0/xdata[i]);
    }
  }
};

// Entropy Problem
template <typename Scalar, uint m, uint n>
struct EntropyProblem : public Problem<Scalar, m, n> {
  EntropyProblem(const EntropyObjective<Scalar, m, n>& obj,
                 const NoConstraints<Scalar, m, n>& nc,
                 const Matrix<Scalar, m, n>& Aeq,
                 const Matrix<Scalar, m, 1>& beq,
                 const GradientFunctor<Scalar, n>& g,
                 const HessianFunctor<Scalar, n>& h) :
      Problem<Scalar, m, n>(obj, nc, Aeq, beq, g, h) { }
};
}  // namespace

TEST(Newton, Quadratic_Programs) {
  const uint n = 3;  // Number of unknowns
  const uint m = 1;  // Number of constraints
  Matrix<double, n, n> Q;
  Matrix<double, n, 1> P;
  double c = 15.0;
  Matrix<double, n, 1> x;

  Q << 5, 0, 0,
       0, 4, 0,
       0, 0, 1.5;
  P << 0.01, 0.5, -0.1;
  x << 100, 100, 100;

  QuadraticObjective<double, n> qobj(Q, P, c);
  const NoConstraints<double, m, n>& nc =
      NoConstraints<double, m, n>::getInstance();
  Matrix<double, m, n> Aeq;
  Matrix<double, m, 1> beq;
  QGradient<double, n> qgradient(Q, P);
  QHessian<double, n> qhessian(Q);

  QuadraticProgram<double, m, n>
      qprogram(qobj, nc, Aeq, beq, qgradient, qhessian);

  Newton<double, m, n> newton(Newton<double, m, n>::UNCONSTRAINED);
  double min_val;
  auto res = newton(qprogram, &x, &min_val);
  LOG(INFO) << "Min_val: " << min_val
            << " Minimizer: [" << x.transpose() << "]"
            << " Result: " << res 
            << std::endl;
  ASSERT_LT(min_val - c, 0.1);
}

TEST(Newton, EqConstrained_Quadratic_Programs) {
  const uint n = 2;  // Number of unknowns
  const uint m = 1;  // Number of constraints
  Matrix<double, n, n> Q;
  Matrix<double, n, 1> P;
  double c = 0.0;
  Matrix<double, n, 1> x;

  Q << 6, 0, 0, 4;
  P << 0.0, 0.0;
  x << 100, 100;

  QuadraticObjective<double, n> qobj(Q, P, c);
  const NoConstraints<double, m, n>& nc =
      NoConstraints<double, m, n>::getInstance();
  Matrix<double, m, n> Aeq;
  Matrix<double, m, 1> beq;
  QGradient<double, n> qgradient(Q, P);
  QHessian<double, n> qhessian(Q);

  Aeq << 1, -1;
  beq << 0;

  QuadraticProgram<double, m, n>
      qprogram(qobj, nc, Aeq, beq, qgradient, qhessian);

  Newton<double, m, n> newton(Newton<double, m, n>::EQUALITY_CONSTRAINED);
  double min_val;
  auto res = newton(qprogram, &x, &min_val);
  LOG(INFO) << "Min_val: " << min_val
            << " Minimizer: [" << x.transpose() << "]"
            << " Result: " << res 
            << std::endl;
  ASSERT_LT(min_val - c, 0.1);
}

TEST(Newton, Unconstrained_Random_Quadratic_Program) {
  const uint n = 16;  // Number of unkowns
  const uint m = 1;  // Unconstrained
  Matrix<double, n, n> Q = Matrix<double, n, n>::Random();
  Matrix<double, n, 1> P = Matrix<double, n, 1>::Zero();
  double c = 5.0;
  Matrix<double, n, 1> x = Matrix<double, n, 1>::Random();

  Q = Q.transpose()*Q;

  QuadraticObjective<double, n> qobj(Q, P, c);
  const NoConstraints<double, m, n>& nc =
      NoConstraints<double, m, n>::getInstance();
  Matrix<double, m, n> Aeq;
  Matrix<double, m, 1> beq;
  QGradient<double, n> qgradient(Q, P);
  QHessian<double, n> qhessian(Q);

  QuadraticProgram<double, m, n>
      qprogram(qobj, nc, Aeq, beq, qgradient, qhessian);

  Newton<double, m, n> newton(Newton<double, m, n>::UNCONSTRAINED);
  double min_val;
  auto res = newton(qprogram, &x, &min_val);
  LOG(INFO) << "Min_val: " << min_val
            << " Minimizer: [" << x.transpose() << "]"
            << " Result: " << res 
            << std::endl;
  ASSERT_LT(min_val - c, 0.1);
}

TEST(Newton, EqConstrained_Random_Quadratic_Programs) {
  const uint n = 16;
  const uint m = 1;
  Matrix<double, n, n> Q =  Matrix<double, n, n>::Random();
  Matrix<double, n, 1> P = Matrix<double, n, 1>::Zero();
  double c = 5.0;
  Matrix<double, n, 1> x;
  x = Matrix<double, n, 1>::Zero();
  x(0, 0) = 1;
  x(n - 1, 0) = -1;

  Q = Q.transpose()*Q;

  QuadraticObjective<double, n> qobj(Q, P, c);
  const NoConstraints<double, m, n>& nc =
      NoConstraints<double, m, n>::getInstance();
  Matrix<double, m, n> Aeq;
  Matrix<double, m, 1> beq;
  QGradient<double, n> qgradient(Q, P);
  QHessian<double, n> qhessian(Q);

  Aeq.setOnes();
  beq = Matrix<double, m, 1>::Zero();

  QuadraticProgram<double, m, n>
      qprogram(qobj, nc, Aeq, beq, qgradient, qhessian);

  Newton<double, m, n> newton(Newton<double, m, n>::EQUALITY_CONSTRAINED);
  double min_val;
  auto res = newton(qprogram, &x, &min_val);
  double rprimal = (Aeq*x-beq).norm();
  LOG(INFO) << "Min_val: " << min_val
            << " Minimizer: [" << x.transpose() << "]"
            << " Result: " << res
            << " Primal Residual: " << rprimal 
            << std::endl;
  ASSERT_LT(min_val - c, 0.1);
  ASSERT_LT(x.sum(), 0.1);
  ASSERT_LT(rprimal, 0.1);
}

TEST(Newton, EqualityConstrained_Entropy) {
  const uint n = 4;
  const uint m = 2;
  double h = 1e-6;
  EntropyObjective<double, m, n> objective;
  const NoConstraints<double, m, n>& nc =
      NoConstraints<double, m, n>::getInstance();
  EntropyHessian<double, n> hessian;
  EntropyGradient<double, n> gradient;
  // NumericalHessian<double, n> hessian(objective, h);
  // SecantGradientFunctor<double, n> gradient(objective, h);
  Matrix<double, m, n> Aeq;
  std::random_device rng;
  std::uniform_int_distribution<> randi(1, n);
  Matrix<double, n, 1> x;

  double* Adata = Aeq.data();
  double* xdata = x.data();
  for (unsigned int i = 0; i < x.size(); i++) xdata[i] = randi(rng);
  uint rank;
  do {
    for (unsigned int i = 0; i < Aeq.size(); i++) Adata[i] = randi(rng);
    Eigen::FullPivLU<Matrix<double, m, n> > lu(Aeq);
    rank = lu.rank();
  } while (rank != m);

  Matrix<double, m, 1> beq = Aeq*x;
  EntropyProblem<double, m, n> problem(objective, nc, Aeq, beq,
                                       gradient, hessian);
  Newton<double, m, n>
      newton(Newton<double, m, n>::INFEASIBLE_EQUALITY_CONSTRAINED);
  newton.options.epsilon_ = 1e-9;
  double min_val;
  auto res = newton(problem, &x, &min_val);
  double rprimal = (Aeq*x-beq).norm();
  LOG(INFO) << "Min_val: " << min_val
            << " Minimizer: [" << x.transpose() << "]"
            << " Result: " << res
            << " Primal Residual: " << rprimal
            << std::endl;
}
}  // solvers
}  // optimo
