// Copyright (C) 2015 The Regents of the University of California (Regents).
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
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include <Eigen/Core>
#include <Eigen/LU>
#include "gtest/gtest.h"

#include "theia/math/l1_solver.h"
#include "theia/util/random.h"

namespace theia {

// A rigged L1 minimization problem with a known output.
//
// [  1  2  3 ]         [ 1 ]
// [  4  5  5 ] * [x] = [ 2 ]
// [  7  8  9 ]         [ 3 ]
// [ 10 11 23 ]         [ 4 ]
//
// For this problem, the L1 minimization should result in
// x = [-0.75, 1.0, 0.0]^t.
TEST(L1Solver, SmallProblem) {
  static const double kTolerance = 1e-8;

  Eigen::MatrixXd lhs(4, 3);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      lhs(i, j) = i * 4 + j;
    }
  }
  lhs(1,2) = 5;

  Eigen::VectorXd rhs(4);
  for (int i = 0; i < 4; i++) {
    rhs(i) = i + 1;
  }

  Eigen::VectorXd solution(3);
  solution.setZero();

  // Recover the code word.
  L1Solver<Eigen::MatrixXd>::Options options;
  options.max_num_iterations = 100;
  L1Solver<Eigen::MatrixXd> l1_solver(options, lhs);
  l1_solver.Solve(rhs, &solution);

  // Verify the solution is near (-0.75, 1.0, 0.0).
  const Eigen::Vector3d gt_solution(-0.75, 1.0, 0.0);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(solution(i), gt_solution(i), kTolerance);
  }

  // Check that the solution is near zero.
  const Eigen::VectorXd residual = lhs * solution - rhs;
  for (int i = 0; i < residual.size(); i++) {
    EXPECT_NEAR(residual(i), 0.0, kTolerance);
  }
}

// This example is taken from the L1-magic library. It is formulated as a
// codeword recovery problem.
TEST(L1Solver, Decoding) {
  RandomNumberGenerator rng(94);
  static const double kTolerance = 1e-8;

  static const int source_length = 256;
  static const int codeword_length = 4 * source_length;
  static const int num_pertubations = 0.2 * codeword_length;

  Eigen::MatrixXd mat(codeword_length, source_length);
  rng.SetRandom(&mat);
  Eigen::VectorXd source_word(source_length);
  rng.SetRandom(&source_word);
  const Eigen::VectorXd code_word = mat * source_word;

  // Apply random pertubations to the initial guess.
  Eigen::VectorXd observation = code_word;
  Eigen::VectorXd pertubations(num_pertubations);
  rng.SetRandom(&pertubations);
  pertubations *= 0.5;
  for (int i = 0; i < num_pertubations; i++) {
    const int rand_entry = rng.RandInt(0, observation.size() - 1);
    observation(rand_entry) = pertubations(i);
  }

  // Recover the code word.
  L1Solver<Eigen::MatrixXd>::Options options;
  options.absolute_tolerance = 1e-8;
  L1Solver<Eigen::MatrixXd> l1_solver(options, mat);
  // Set the initial noisy guess.
  Eigen::VectorXd solution =
      (mat.transpose() * mat).inverse() * mat.transpose() * observation;
  l1_solver.Solve(observation, &solution);

  // Verify the solution.
  const Eigen::VectorXd residual = mat * solution - code_word;
  for (int i = 0; i < residual.size(); i++) {
    EXPECT_NEAR(residual(i), 0.0, kTolerance);
  }

}

}  // namespace theia
