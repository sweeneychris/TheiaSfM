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

#include <glog/logging.h>
#include "gtest/gtest.h"
#include "optimo/core/objects.h"
#include "optimo/core/numerical_gradient.h"

namespace optimo {
using Eigen::Matrix;

namespace objects {
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
    return sum(0, 0) + c;
  }

 protected:
  Matrix<Scalar, n, n> Q;
  Matrix<Scalar, n, 1> P;
  Scalar c;
};
}  // namespace

TEST(Numerical_Gradient, SecantGradientFunctor_Computation) {
  const uint n = 3;
  const float h = 1e-6;
  Matrix<float, n, n> Q;
  Matrix<float, n, 1> P;
  Matrix<float, n, 1> x;
  float c = 5.0;
  Q << 2, 0, 0,
       0, 5, 0,
       0, 0, 8;
  P << 1, 1, 1;
  x << 1, 0, 0;
  QuadraticFunctor<float, n> quadratic(Q, P, c);
  float r = quadratic(x);
  SecantGradientFunctor<float, 3> gradient(quadratic, h);

  // Case 1
  Matrix<float, n, 1> g_gt = 2*Q*x + P;
  Matrix<float, n, 1> g;
  gradient(x, &g);
  Matrix<float, n, 1> error = g_gt - g;
  LOG(INFO) << "Comp. Gradient: [" << g.transpose() << "]'";
  LOG(INFO) << "Gradient GT: [" << g_gt.transpose() << "]'";
  LOG(INFO) << "Vector Error: [" << error.transpose() << "]'";
  ASSERT_LT(error.norm(), 1.0f);

  // Case 2 x = 0
  x << 0, 0, 0;
  g_gt = P;
  gradient(x, &g);
  error = g_gt - g;
  LOG(INFO) << "Comp. Gradient: [" << g.transpose() << "]'";
  LOG(INFO) << "Gradient GT: [" << g_gt.transpose() << "]'";
  LOG(INFO) << "Vector Error: [" << error.transpose() << "]'";
  ASSERT_LT(error.norm(), 1.0f);

  // Double PRecision
  Matrix<double, n, n> Qd;
  Matrix<double, n, 1> Pd;
  Matrix<double, n, 1> xd;
  Qd << 2, 0, 0,
        0, 5, 0,
        0, 0, 8;
  Pd << 1, 1, 1;
  xd << 1, 0, 0;
  QuadraticFunctor<double, n> quadraticd(Qd, Pd, static_cast<double>(c));
  SecantGradientFunctor<double, 3> gradientd(quadraticd,
                                             static_cast<double>(h));

  // Case 1
  Matrix<double, n, 1> g_gtd = 2*Qd*xd + Pd;
  Matrix<double, n, 1> gd;
  gradientd(xd, &gd);
  Matrix<double, n, 1> errord = g_gtd - gd;
  LOG(INFO) << "Comp. Gradient: [" << gd.transpose() << "]'";
  LOG(INFO) << "Gradient GT: [" << g_gtd.transpose() << "]'";
  LOG(INFO) << "Vector Error: [" << errord.transpose() << "]'";
  ASSERT_LT(errord.norm(), 1.0f);

  // Case 2 x = 0
  xd << 0, 0, 0;
  g_gtd = Pd;
  gradientd(xd, &gd);
  errord = g_gtd - gd;
  LOG(INFO) << "Comp. Gradient: [" << gd.transpose() << "]'";
  LOG(INFO) << "Gradient GT: [" << g_gtd.transpose() << "]'";
  LOG(INFO) << "Vector Error: [" << errord.transpose() << "]'";
  ASSERT_LT(errord.norm(), 1.0f);
}
// TODO(vfragoso): Add tests for numerical gradient LS
}  // objects
}  // optimo
