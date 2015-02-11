// Copyright (C) 2013  Victor Fragoso <vfragoso@cs.ucsb.edu>
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

#include "gtest/gtest.h"
#include "optimo/core/objects.h"
#include "optimo/core/numerical_gradient.h"
#include "optimo/core/numerical_hessian.h"
#include <vector>

namespace optimo {

using Eigen::Matrix;
using std::vector;

namespace {
template <typename Scalar, uint n>
class QuadraticFunctor : public ObjectiveFunctor<Scalar, n> {
 public:
  QuadraticFunctor(const Matrix<Scalar, n, n>& _Q,
                   const Matrix<Scalar, n, 1>& _P,
                   const Scalar& _c) {
    Q = _Q;
    P = _P;
    c = _c;
  }

  ~QuadraticFunctor(void) { }
  
  virtual Scalar operator()(const Matrix<Scalar, n, 1>& x) const {
    Matrix<Scalar, 1, 1> sum = x.transpose()*Q*x + P.transpose()*x;
    return sum(0,0) + c;
  }
  
 protected:
  Matrix<Scalar, n, n> Q;
  Matrix<Scalar, n, 1> P;
  Scalar c;
};

template <typename Scalar, uint n>
class LinearConstraint : public ConstraintFunctor<Scalar, n> {
 public:
  LinearConstraint(const Matrix<Scalar, 1, n>& _A) {
    A = _A;
  }

  ~LinearConstraint(void) { }

  virtual Scalar operator()(const Matrix<Scalar, n, 1>& x) const {
    Matrix<Scalar, 1, 1> t = A*x;
    return static_cast<Scalar>(t(0,0));
  }

 protected:
  Matrix<Scalar, 1, n> A;
};

template <typename Scalar, uint m, uint n>
class VectorConstraints : public Constraints<Scalar, m, n> {
 public:
  VectorConstraints(const vector<ConstraintFunctor<Scalar, n>* >& _vec) {
    vec = _vec;
  }

  virtual ~VectorConstraints(void) { }

  virtual ConstraintFunctor<Scalar, n>&
  operator[](const unsigned int i) {
    return *vec[i];
  }
  
 protected:
  vector<ConstraintFunctor<Scalar, n>* > vec;
};

template <typename Scalar, uint m, uint n>
struct DummyProblem : public Problem<Scalar, m, n> {
  DummyProblem(const ObjectiveFunctor<Scalar, n>& obj,
               const GradientFunctor<Scalar, n>& g,
               const HessianFunctor<Scalar, n>& h) :
      Problem<Scalar, m, n>(obj, g, h) { }

  virtual ~DummyProblem(void) { }
};
}  // namespace

TEST(Objects, ObjectiveFunctor_Creation) {
  const uint n = 3;
  Matrix<float, n, n> Q;
  Matrix<float, n, 1> P;
  Matrix<float, n, 1> x;
  float c = 5.0;
  Q << 1, 0, 0, 0 , 1, 0, 0, 0, 1;
  P << 1, 1, 1;
  x << 1, 0, 0;
  QuadraticFunctor<float, n> quadratic(Q, P, c);
  float r = quadratic(x);
  float expected = 7;
  ASSERT_EQ(r, expected);
}

TEST(Objects, ConstraintFunctor_Creation) {
  const uint n = 3;
  Matrix<float, n, 1> x;
  Matrix<float, 1, n> A;
  A << 1, 1, 1;
  x << 1, 0, 0;
  LinearConstraint<float, n> constraint(A);
  float r = constraint(x);
  float expected = 1;
  ASSERT_EQ(r, expected);
}

TEST(Objects, Constraints_Creation) {
  const uint n = 3;
  const uint m = 2;
  Matrix<float, n, 1> x;
  Matrix<float, 1, n> A1, A2;
  A1 << 1, 1, 1;
  A2 << 0, 0, 0;
  x << 1, 0, 0;
  LinearConstraint<float, n> constraint1(A1);
  LinearConstraint<float, n> constraint2(A2);
  vector<ConstraintFunctor<float, n>*> vec_constraints;
  vec_constraints.push_back(&constraint1);
  vec_constraints.push_back(&constraint2);
  VectorConstraints<float, m, n> constraints(vec_constraints);

  float r1 = constraints[0](x);
  float r2 = constraints[1](x);

  float expected1 = 1;
  float expected2 = 0;

  ASSERT_EQ(r1, expected1);
  ASSERT_EQ(r2, expected2);
}

TEST(Objects, Problem_Creation) {
  const uint m = 0;
  const uint n = 3;
  double h = 1e-6;
  Matrix<double, n, n> Q;
  Matrix<double, n, 1> P;
  Matrix<double, n, 1> x;
  double c = 5.0;
  Q << 1, 0, 0, 0 , 1, 0, 0, 0, 1;
  P << 1, 1, 1;
  x << 1, 0, 0;

  QuadraticFunctor<double, n> quadratic(Q, P, c);
  optimo::SecantGradientFunctor<double, n> qgradient(quadratic, h);
  optimo::NumericalHessian<double, n> qhessian(quadratic, h);
  
  DummyProblem<double, m, n> problem(quadratic, qgradient, qhessian);
  // TODO(vfragoso): Add assertions here!
}
}  // optimo

