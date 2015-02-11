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
#include "optimo/core/objects_ls.h"
#include "optimo/solvers/gradient_descent.h"
#include "optimo/core/numerical_gradient.h"
#include "optimo/core/numerical_hessian.h"
#include <glog/logging.h>
#include "gtest/gtest.h"
#include <Eigen/LU>

namespace optimo {
namespace solvers {
namespace {
using Eigen::Matrix;
using Eigen::Dynamic;
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


// Quadratic Objective LS
template <typename Scalar>
class QuadraticObjectiveLS : public ObjectiveFunctorLS<Scalar> {
 public:
  QuadraticObjectiveLS(const Matrix<Scalar, Dynamic, Dynamic>& _Q,
                       const Matrix<Scalar, Dynamic, 1>& _P,
                       const Scalar& _c) {
    Q = _Q;
    P = _P;
    c = _c;
  }

  ~QuadraticObjectiveLS(void) { }
  
  virtual Scalar operator()(const Matrix<Scalar, Dynamic, 1>& x) const {
    Matrix<Scalar, 1, 1> sum = x.transpose()*Q*x + P.transpose()*x;
    return sum(0,0) + c;
  }
  
 protected:
  Matrix<Scalar, Dynamic, Dynamic> Q;
  Matrix<Scalar, Dynamic, 1> P;
  Scalar c;
};

// Quadratic Program
template <typename Scalar, uint m, uint n>
class QuadraticProgram : public Problem<Scalar, m, n> {
 public:
  QuadraticProgram(const QuadraticObjective<Scalar, n>& qobj,
                   const GradientFunctor<Scalar, n>& g,
                   const HessianFunctor<Scalar, n>& h) :
      Problem<Scalar, m, n>(qobj, g, h) {
  }

  virtual ~QuadraticProgram(void) { }
};

// Quadratic Program LS
template <typename Scalar>
class QuadraticProgramLS : public ProblemLS<Scalar> {
 public:
  QuadraticProgramLS(const QuadraticObjectiveLS<Scalar>& qobj,
                   const GradientFunctorLS<Scalar>& g,
                   const HessianFunctorLS<Scalar>& h) :
      ProblemLS<Scalar>(qobj, g, h) {
  }

  virtual ~QuadraticProgramLS(void) { }
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

// Gradient Functor LS
template <typename Scalar>
class QGradientLS : public GradientFunctorLS<Scalar> {
 public:
  QGradientLS(const Matrix<Scalar, Dynamic, Dynamic>& _Q,
              const Matrix<Scalar, Dynamic, 1>& _P) {
    Q = _Q;
    P = _P;
  }

  virtual void operator()(const Matrix<Scalar, Dynamic, 1>& x,
                          Matrix<Scalar, Dynamic, 1>* g) const {
    *g = 2*Q*x + P;
  }

 protected:
  Matrix<Scalar, Dynamic, Dynamic> Q;
  Matrix<Scalar, Dynamic, 1> P;
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

// Hessian Functor
template <typename Scalar>
class QHessianLS : public HessianFunctorLS<Scalar> {
 public:
  QHessianLS(const Matrix<Scalar, Dynamic, Dynamic>& _Q) {
    Q = _Q;
  }

  virtual void operator()(const Matrix<Scalar, Dynamic, 1>& x,
                          Matrix<Scalar, Dynamic, Dynamic>* h) const {
    *h = 2*Q;
  }
  
 protected:
  Matrix<Scalar, Dynamic, Dynamic> Q;
};
}  // namespace

TEST(GradientDescent, Quadratic_Programs) {
  const uint n = 3;
  const uint m = 0;
  Matrix<double, n, n> Q;
  Matrix<double, n, 1> P;
  double c = 15.0;
  Matrix<double, n, 1> x;

  Q << 5, 0, 0,
       0, 4, 0,
       0, 0, 1.5;
  P << 0.01, 0.5, -0.1;
  x << 101, 101, 101;

  QuadraticObjective<double, n> qobj(Q, P, c);
  QGradient<double, n> qgradient(Q, P);
  QHessian<double, n> qhessian(Q);

  QuadraticProgram<double, m, n>
      qprogram(qobj, qgradient, qhessian);

  GradientDescent<double, m, n> gdescent;
  double min_val;
  auto res = gdescent(qprogram, &x, &min_val);

  LOG(INFO) << "Min_val: " << min_val
            << " Minimizer: [" << x.transpose() << "]"
            << " Result: " << res 
            << std::endl;
  ASSERT_LT(min_val - c, 0.1);
}

TEST(GradientDescent, Quatratic_ProgramLS) {
  const uint n = 3;
  const uint m = 0;
  Matrix<double, Dynamic, Dynamic> Q(n, n);
  Matrix<double, Dynamic, 1> P(n);
  double c = 15.0;
  Matrix<double, Dynamic, 1> x(n);

  Q << 5, 0, 0,
       0, 4, 0,
       0, 0, 1.5;
  P << 0.01, 0.5, -0.1;
  x << 101, 101, 101;

  QuadraticObjectiveLS<double> qobj(Q, P, c);
  QGradientLS<double> qgradient(Q, P);
  QHessianLS<double> qhessian(Q);

  QuadraticProgramLS<double>
      qprogram(qobj, qgradient, qhessian);

  GradientDescentLS<double> gdescent(n);
  double min_val;
  auto res = gdescent(qprogram, &x, &min_val);

  LOG(INFO) << "Min_val: " << min_val
            << " Minimizer: [" << x.transpose() << "]"
            << " Result: " << res 
            << std::endl;
  ASSERT_LT(min_val - c, 0.1);
}
}  // solvers
}  // optimo
