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

#ifdef WITH_SUITESPARSE
#include <glog/logging.h>
#include "gtest/gtest.h"
#include "optimo/solvers/sparse_primal_dual_qp.h"

namespace optimo {
namespace solvers {
using Eigen::Matrix;
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::Vector2d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Dynamic;

TEST(Sparse_Primal_Dual_Interior_Point_QP, NoData) {
  const uint n = 3;
  const uint m = 4;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;
  SparsePrimalDualQP<double>::Params params(n, m, p);
  params.d(n - 1) = 1.0;
  params.beq(0) = 1.0;
  SparsePrimalDualQP<double> qp_solver;
  Matrix<double, Dynamic, 1> y(l);
  y.setConstant(1.0);
  y(0) = 0.0;
  y(1) = 2.0;
  double min_value;
  TERMINATION_TYPE r = qp_solver(&params, &y, &min_value);
  ASSERT_EQ(r, INVALID_ARGUMENTS);
}

TEST(Sparse_Primal_Dual_Interior_Point_QP, InvalidArguments) {
  const uint n = 3;
  const uint m = 4;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;
  SparsePrimalDualQP<double>::Params params(n, m, p);
  SparsePrimalDualQP<double> qp_solver;
  qp_solver.options.alpha_ = -1.0;
  qp_solver.options.beta_ = -1.0;
  Matrix<double, Dynamic, 1> y(l);
  double min_value;
  TERMINATION_TYPE r = qp_solver(&params, &y, &min_value);
  ASSERT_EQ(r, INVALID_ARGUMENTS);
}

TEST(Sparse_Primal_Dual_Interior_Point_QP, Small_QP) {
  const uint n = 2;
  const uint m = 5;  // To ensure x >= 0
  const uint p = 1;  // Dummy eq. constraint
  const uint l = n + m + p;
  SparsePrimalDualQP<double>::Params params(n, m, p);
  // To allocate memory for F (optional). This speeds things up.
  params.reserve();
  // Q = [ 1 -1;
  //      -1 2]
  // d = [-2; -6]
  // Ain = [1 1;
  //       -1 2;   <= b = [2; 2; 3]
  //        2 1]*x 
  // x >= 0
  // Aeq = [1 1] beq = [1]

  // Setting Q
  params.setQElement(0, 0, 1.0);
  params.setQElement(0, 1, -1.0);
  params.setQElement(1, 0, -1.0);
  params.setQElement(1, 1, 2.0);

  // Setting Ain
  params.setAinElement(0, 0, 1.0);
  params.setAinElement(0, 1, 1.0);
  params.setAinElement(1, 0, -1.0);
  params.setAinElement(1, 1, 2.0);
  params.setAinElement(2, 0, 2.0);
  params.setAinElement(2, 1, 1.0);
  params.setAinElement(3, 0, -1.0);
  params.setAinElement(4, 1, -1.0);

  // Setting Aeq
  params.setAeqElement(0, 0, 1.0);
  params.setAeqElement(0, 1, 1.0);

  // Setting d
  Matrix<double, Dynamic, 1>& d = params.d;
  d(0) = -2;
  d(1) = -6;

  // Setting beq
  Matrix<double, Dynamic, 1>& beq = params.beq;
  beq(0) = 1.0;

  // Setting bin
  Matrix<double, Dynamic, 1>& bin = params.bin;
  bin.setConstant(0.0);
  bin(0) = 2.0;
  bin(1) = 2.0;
  bin(2) = 3.0;

  // Initial point
  Matrix<double, Dynamic, 1> y(l);
  y(0) = 0.5;
  y(1) = 0.5;
  
  SparsePrimalDualQP<double> qp_solver;
  double min_value;

  auto res = qp_solver(&params, &y, &min_value);
  ASSERT_EQ(res, 0);

  VLOG(1) << y.transpose();
  ASSERT_NEAR(y(0), 0.0, 1e-6);
  ASSERT_NEAR(y(1), 1.0, 1e-6);
}

TEST(Sparse_Primal_Dual_Interior_Point_QP, SingleClass_SVM) {
  // Large PSD
  const int dim = 250;
  Matrix<double, Dynamic, Dynamic> A(dim, dim);
  A.setRandom();
  Matrix<double, Dynamic, Dynamic> Q = A.transpose()*A;
  const int n = dim;
  const int m = 2*n;
  const int p = 1;
  const uint l = n + m + p;
  const double nu = 0.5;
  const double c = 1.0/(nu*dim);
  SparsePrimalDualQP<double>::Params params(n, m, p);
  // To allocate memory for F (optional). This speeds things up.
  params.reserve();
  SparsePrimalDualQP<double> qp_solver;
  Matrix<double, Dynamic, 1>& d = params.d;
  Matrix<double, Dynamic, 1>& beq = params.beq;
  Matrix<double, Dynamic, 1>& bin = params.bin;
  Matrix<double, Dynamic, 1> y(l);

  d.setConstant(0.0);
  beq(0) = 1.0;
  bin.setConstant(0.0);
  for (int i = 0; i < dim; i++) {
    for (int j = i; j < dim; j++) {
      params.setQElement(i, j, Q(i, j));
      params.setQElement(j, i, Q(j, i));
    }
    params.setAeqElement(0, i, 1.0);
    params.setAinElement(i, i, -1.0);
    params.setAinElement(n + i, i, 1.0);
    bin(n + i) = c;
  }
  y.block(0, 0, n, 1).setConstant(1.0/dim);

  double min_value;
  auto res = qp_solver(&params, &y, &min_value);
  ASSERT_EQ(res, SOLVED);
  VLOG(1) << y.block(0, 0, n, 1).transpose();
}
}  // solvers
}  // optimo
#endif  // WITH_SUITESPARSE
