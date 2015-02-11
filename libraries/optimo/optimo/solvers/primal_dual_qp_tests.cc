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

#include "optimo/solvers/primal_dual_qp.h"
#include <glog/logging.h>
#include "gtest/gtest.h"

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::Dynamic;

TEST(Primal_Dual_QP, Constrained_LeastSquares) {
  // Q = [ 1 -1;
  //      -1 2]
  // d = [-2; -6]
  // Ain = [1 1;
  //       -1 2;   <= b = [2; 2; 3]
  //        2 1]*x 
  // x >= 0
  const uint n = 2;  // n number of unknowns
  const uint m = 5;  // m number of inequalities
  const uint p = 1;  // p number of equalities
  PrimalDualQP<double, n, m, p> qp_solver;
  qp_solver.options.max_iter_ = 100;
  PrimalDualQP<double, n, m, p>::Params params;
  
  // Setting matrix Q
  auto& Q = params.Q;
  Q *= 0.0;  // Setting to zero
  Q(0, 0) = 1; Q(0, 1) = -1;
  Q(1, 0) = -1; Q(1, 1) = 2;
  VLOG(1) << "Q: \n" << Q;

  // Setting Ain
  auto& Ain = params.Ain;
  Ain.setConstant(0);
  Ain(0, 0) = 1; Ain(0, 1) = 1;
  Ain(1, 0) = -1; Ain(1, 1) = 2;
  Ain(2, 0) = 2; Ain(2, 1) = 1;
  Ain(3, 0) = -1;  // for x >= 0
  Ain(4, 1) = -1;  // for x >= 0
  VLOG(1) << "Ain: \n" << Ain;

  // Setting Aeq
  auto& Aeq = params.Aeq;
  Aeq(0, 0) = 1; Aeq(0, 1) = 1;
  VLOG(1) << "Aeq: \n" << Aeq;

  // Setting bin
  auto& bin = params.bin;
  bin.setConstant(0);
  bin(0) = 2;
  bin(1) = 2;
  bin(2) = 3;
  VLOG(1) << "bin: \n" << bin;

  auto& d = params.d;
  d(0) = -2;
  d(1) = -6;
  VLOG(1) << "d: \n" << d;

  auto& beq = params.beq;
  beq(0) = 1;
  VLOG(1) << "beq: \n" << beq;
  
  Matrix<double, n, 1> x;
  x(0) = 0.5;
  x(1) = 0.5;
  double min_value;
  auto res = qp_solver(params, &x, &min_value);
  ASSERT_EQ(res, 0);
  VLOG(1) << "Minimizer: " << x.transpose()
          << " min_value=" << min_value;

  // Values calculated w/ MATLAB
  ASSERT_NEAR(x(0), 0.0, 1e-3);
  ASSERT_NEAR(x(1), 1.0, 1e-3);
  ASSERT_NEAR(min_value, -5.0, 1e-3);
}

TEST(Primal_Dual_QP, Constrained_LeastSquares2) {
  //       Q               d
  // 1090.27 00430.3     -770.883
  // 00430.3 261.802     -354.634
  // Ain       bin
  // -1 00      0
  // 00 -1      0
  // 01 00 <=   0.5
  // 00 01      1
  //
  // Aeq       beq
  // 1 1    =   1
  const uint n = 2;  // n number of unknowns
  const uint m = 4;  // m number of inequalities
  const uint p = 1;  // p number of equalities
  PrimalDualQP<double, n, m, p> qp_solver;
  qp_solver.options.max_iter_ = 100;
  PrimalDualQP<double, n, m, p>::Params params;
  
  // Setting matrix Q
  auto& Q = params.Q;
  Q *= 0.0;  // Setting to zero
  Q(0, 0) = 1090.27; Q(0, 1) = 430.3;
  Q(1, 0) = 430.3; Q(1, 1) = 261.802;
  VLOG(1) << "Q: \n" << Q;

  auto& d = params.d;
  d(0) = -770.883;
  d(1) = -354.634;
  VLOG(1) << "d: " << d.transpose();

  // Setting Ain
  auto& Ain = params.Ain;
  Ain.setConstant(0);
  Ain(0, 0) = -1; 
  Ain(1, 1) = -1;
  Ain(2, 0) = 1;
  Ain(3, 1) = 1;
  VLOG(1) << "Ain: \n" << Ain;

  // Setting Aeq
  auto& Aeq = params.Aeq;
  Aeq(0, 0) = 1;
  Aeq(0, 1) = 1;
  VLOG(1) << "Aeq: \n" << Aeq;

  // Setting bin
  auto& bin = params.bin;
  bin.setConstant(0);
  bin(2) = 0.5;
  bin(3) = 1.0;
  VLOG(1) << "bin: \n" << bin.transpose();

  auto& beq = params.beq;
  beq(0) = 1;
  VLOG(1) << "beq: " << beq.transpose();
  
  Matrix<double, n, 1> x;
  x(0) = 0.5;
  x(1) = 0.5;
  double min_value;
  auto res = qp_solver(params, &x, &min_value);

  VLOG(1) << "Minimizer: " << x.transpose()
          << " min_value=" << min_value;

  // Values calculated w/ MATLAB
  ASSERT_EQ(res, 0);
  ASSERT_NEAR(x(0), 0.5, 1e-3);
  ASSERT_NEAR(x(1), 0.5, 1e-3);
}

TEST(Primal_Dual_QP, SingleClass_SVM) {
  const int dim = 250;
  const int n = dim;
  const int m = 2*n;
  const int p = 1;
  const int l = n + m + p;
  const double nu = 0.5;
  const double c = 1.0/(nu*dim);
  PrimalDualQP<double, n, m, p> qp_solver;
  PrimalDualQP<double, n, m, p>::Params params;
  Matrix<double, n, 1> x;
  Matrix<double, Dynamic, Dynamic> A(dim, dim);
  A.setRandom();
  params.Q = A.transpose()*A;
  params.d.setConstant(0.0);
  params.Aeq.setConstant(1.0);
  params.beq(0) = 1.0;
  params.bin.block(0, 0, n, 1).setConstant(0.0);
  params.bin.block(n, 0, n, 1).setConstant(c);
  x.setConstant(1.0/dim);
  for (int i = 0; i < dim; i++) {
    params.Ain(i, i) = -1.0;
    params.Ain(n + i, i) = 1.0;
  }
  double min_value;
  auto res = qp_solver(params, &x, &min_value);
  ASSERT_EQ(res, SOLVED);
  VLOG(1) << x.transpose();
}
}  // solvers
}  // optimo
