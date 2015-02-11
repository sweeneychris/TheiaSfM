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

#include "optimo/solvers/primal_dual_lp.h"
#include <glog/logging.h>
#include "gtest/gtest.h"

namespace optimo {
namespace solvers {

using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Dynamic;

TEST(Primal_Dual_Interior_Point_LP, NoData) {
  const uint n = 3;
  const uint m = 4;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;
  PrimalDualLP<double, n, m, p> lp_solver;
  PrimalDualLP<double, n, m, p>::Params params;
  Matrix<double, Dynamic, 1> y(l);
  // y = [x lambdas nus]'
  y(0) = 0.0;
  y(1) = 2.0;
  Matrix<double, Dynamic, 1>& c = params.c;
  Matrix<double, Dynamic, 1>& beq = params.beq;
  c(n - 1) = 1.0;
  beq(0) = 1.0;
  double min_value;
  auto r = lp_solver(params, &y, &min_value);
  ASSERT_EQ(r, INVALID_ARGUMENTS);
}


TEST(Primal_Dual_Interior_Point_LP, InvalidArguments) {
  const uint n = 3;
  const uint m = 4;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;
  PrimalDualLP<double, n, m, p> lp_solver;
  PrimalDualLP<double, n, m, p>::Params params;
  Matrix<double, Dynamic, 1> y(l);
  lp_solver.options.alpha_ = -1;
  lp_solver.options.beta_ = -1;
  double min_value;
  auto r = lp_solver(params, &y, &min_value);
  ASSERT_EQ(r, INVALID_ARGUMENTS);
}

TEST(Primal_Dual_Interior_Point_LP, Small_Game) {
  const uint n = 3;
  const uint m = 4;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;
  PrimalDualLP<double, n, m, p> lp_solver;
  PrimalDualLP<double, n, m, p>::Params params;
  Matrix<double, Dynamic, 1> y(l);
  y(0) = 0.5;
  y(1) = 0.5;
  y(2) = 1.0/3.0;
  // Simple Game
  // G = [3 0; -1 1];
  // Ain = [-G -ones(size(G,1), 1); -eye(size(G,2)) zeros(size(G,2), 1)];
  // bin = zeros(size(Ain,1), 1);
  // Aeq = [ones(1, size(G,2)) 0];
  // beq = 1;
  // c = [zeros(1, size(G,2)) 1]';
  Matrix<double, Dynamic, 1>& beq = params.beq;
  Matrix<double, Dynamic, 1>& c = params.c;

  c.setConstant(0.0);
  c(2) = 1.0;
  beq(0) = 1.0;

  // Setting Aeq (Prob Simplex)
  params.Aeq(0, 0) = 1.0;
  params.Aeq(0, 1) = 1.0;
  params.Aeq(0, 2) = 0.0;

  // Setting Ain (Inequality)
  params.Ain(0, 0) = -3.0;
  params.Ain(0, 2) = -1.0;
  params.Ain(1, 0) = 1.0;
  params.Ain(1, 1) = -1.0;
  params.Ain(1, 2) = -1.0;
  params.Ain(2, 0) = -1.0;
  params.Ain(3, 1) = -1.0;
  params.bin.setConstant(0.0);

  double min_value;
  auto r = lp_solver(params, &y, &min_value);
  ASSERT_EQ(SOLVED, r);

  Matrix<double, Dynamic, 1> x_star(n);
  x_star(0) = 0.2;
  x_star(1) = 0.8;
  x_star(2) =-0.6;
  double error = (y.block(0, 0, n, 1) - x_star).norm();
  ASSERT_LT(error, 1e-3);
}

TEST(Primal_Dual_Interior_Point_LP, Random_LP) {
  // c = [-0.7818 0.9494 0.0698]';
  // Ain = [-0.7348   -0.0256    1.4016;
  //        -0.8550   -1.3465    1.0384;
  //         1.9888   -0.3194   -0.0228];
  // bin = zeros(m, 1);
  // Aeq = [0.9293    0.3706   -0.9302];
  // beq = 0.8957;

  const uint n = 3;
  const uint m = 3;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;

  PrimalDualLP<double, n, m, p> lp_solver;
  PrimalDualLP<double, n, m, p>::Params params;
  Matrix<double, Dynamic, 1> y(l);
  y(0) = -0.0255;
  y(1) = 1.406;
  y(2) = -1.0048;

  Matrix<double, Dynamic, 1>& bin = params.bin;
  Matrix<double, Dynamic, 1>& beq = params.beq;
  Matrix<double, Dynamic, 1>& c = params.c;

  bin.setConstant(0.0);

  c(0) = -0.7818;
  c(1) =  0.9494;
  c(2) =  0.0698;
  beq(0) = 0.8957;

  // Setting Aeq (Eq. Constraints)
  params.Aeq(0, 0) = 0.9293;
  params.Aeq(0, 1) = 0.3706;
  params.Aeq(0, 2) = -0.9302;

  // Setting Ain (Inequality Constraints)
  params.Ain(0, 0) = -0.7348;
  params.Ain(0, 1) = -0.0256;
  params.Ain(0, 2) = 1.4016;
  params.Ain(1, 0) = -0.8550;
  params.Ain(1, 1) = -1.3465;
  params.Ain(1, 2) = 1.0384;
  params.Ain(2, 0) = 1.9888;
  params.Ain(2, 1) = -0.3194;
  params.Ain(2, 2) = -0.0228;

  double min_value;
  auto r = lp_solver(params, &y, &min_value);
  ASSERT_EQ(SOLVED, r);

  Matrix<double, Dynamic, 1> x_star(n);
  x_star(0) = -0.1967;
  x_star(1) = -1.1104;
  x_star(2) = -1.6018;
  double error = (y.block(0, 0, n, 1)- x_star).norm();
  ASSERT_LT(error, 1e-3);
}

TEST(Primal_Dual_Interior_Point_LP, Game2) {
  // G = [1 3 0;
  //      6 2 7];
  // Ain = [-G -ones(size(G,1), 1); -eye(size(G,2)) zeros(size(G,2), 1)];
  // bin = zeros(size(Ain,1), 1);
  // Aeq = [ones(1, size(G,2)) 0];
  // beq = 1;
  // c = [zeros(1, size(G,2)) 1]';
  const uint n = 4;
  const uint m = 5;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;

  Matrix<double, Dynamic, 1> y(l);
  PrimalDualLP<double, n, m, p> lp_solver;
  PrimalDualLP<double, n, m, p>::Params params;

  y.block(0, 0, n, 1).setConstant(0.25);
  Matrix<double, Dynamic, 1>& beq = params.beq;
  Matrix<double, Dynamic, 1>& c = params.c;

  c(n - 1, 0) = 1.0;
  beq(0, 0) = 1.0;
  params.Aeq.block(0, 0, p, n - 1).setConstant(1.0);

  // Setting Ain
  params.Ain(0, 0) = -1.0;
  params.Ain(0, 1) = -3.0;
  params.Ain(0, 3) = -1.0;
  params.Ain(1, 0) = -6.0;
  params.Ain(1, 1) = -2.0;
  params.Ain(1, 2) = -7.0;
  params.Ain(1, 3) = -1.0;
  params.Ain(2, 0) = -1.0;
  params.Ain(3, 1) = -1.0;
  params.Ain(4, 2) = -1.0;
  params.bin.setConstant(0.0);

  double min_value;
  auto r = lp_solver(params, &y, &min_value);
  ASSERT_EQ(SOLVED, r);

  Matrix<double, Dynamic, 1> x_star(n);
  x_star(0) = 0.1667;
  x_star(1) = 0.8333;
  x_star(2) = 0.0;
  x_star(3) = -2.6667;
  double error = (y.block(0, 0, n, 1) - x_star).norm();
  ASSERT_LT(error, 1e-3);
}

// RPS: Rock Paper and Scissors Game
TEST(Primal_Dual_Interior_Point_LP, RPS_Game) {
  //   G = [0  1 -1;
  //     -1  0  1;
  //      1 -1  0];
  // Ain = [-G -ones(size(G,1), 1); -eye(size(G,2)) zeros(size(G,2), 1)];
  // bin = zeros(size(Ain,1), 1);
  // Aeq = [ones(1, size(G,2)) 0];
  // beq = 1;
  // c = [zeros(1, size(G,2)) 1]';
  const uint n = 4;
  const uint m = 6;  // To ensure x >= 0
  const uint p = 1;
  const uint l = n + m + p;

  PrimalDualLP<double, n, m, p> lp_solver;
  PrimalDualLP<double, n, m, p>::Params params;
  Matrix<double, Dynamic, 1> y(l);
  y.block(0, 0, n, 1).setConstant(0.25);

  Matrix<double, Dynamic, 1>& beq = params.beq;
  Matrix<double, Dynamic, 1>& c = params.c;

  c(n - 1) = 1.0;
  beq(0) = 1.0;
  params.Aeq.block(0, 0, p, n - 1).setConstant(1.0);

  // Setting Ain
  params.Ain(0, 1) = -1.0;
  params.Ain(0, 2) = 1.0;
  params.Ain(0, 3) = -1.0;
  params.Ain(1, 0) =  1.0;
  params.Ain(1, 2) = -1.0;
  params.Ain(1, 3) = -1.0;
  params.Ain(2, 0) = -1.0;
  params.Ain(2, 1) =  1.0;
  params.Ain(2, 3) = -1.0;
  params.Ain(3, 0) = -1.0;
  params.Ain(4, 1) = -1.0;
  params.Ain(5, 2) = -1.0;

  double min_value;
  auto r = lp_solver(params, &y, &min_value);
  ASSERT_EQ(SOLVED, r);

  Matrix<double, Dynamic, 1> x_star(n);
  x_star(0) = 0.3333;
  x_star(1) = 0.3333;
  x_star(2) = 0.3333;
  x_star(3) = 0.0;
  double error = (y.block(0, 0, n, 1) - x_star).norm();
  ASSERT_LT(error, 1e-3);
}
}  // solvers
}  // optimo
