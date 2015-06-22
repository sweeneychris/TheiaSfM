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
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
#include <glog/logging.h>

#include "gtest/gtest.h"
#include "theia/math/matrix/dominant_eigensolver.h"
#include "theia/math/matrix/linear_operator.h"

namespace theia {

TEST(DominantEigensolver, LargestEigenvalue) {
  static const int kNumDimensions = 10;
  static const double kEigenvalueTolerance = 1e-4;
  static const double kEigenvectorTolerance = 1e-4;

  Eigen::MatrixXd matrix(kNumDimensions, kNumDimensions);
  matrix.setRandom();
  matrix = (matrix.transpose() * matrix).eval();
  DenseLinearOperator linear_operator(matrix);
  DominantEigensolver::Options options;
  DominantEigensolver eigensolver(options, linear_operator);

  // Compute with power iterations.
  Eigen::VectorXd eigenvector;
  double eigenvalue;
  EXPECT_TRUE(eigensolver.Compute(&eigenvalue, &eigenvector));

  // Get ground truth.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gt_eigensolver(matrix);
  const Eigen::VectorXd gt_eigenvector =
      gt_eigensolver.eigenvectors().rightCols<1>();
  const double gt_eigenvalue = gt_eigensolver.eigenvalues().reverse()(0);

  EXPECT_NEAR(eigenvalue, gt_eigenvalue, kEigenvalueTolerance);
  EXPECT_LT(1.0 - std::abs(eigenvector.dot(gt_eigenvector)),
            kEigenvectorTolerance);
}

TEST(DominantEigensolver, SmallestEigenvalue) {
  static const int kNumDimensions = 10;
  static const double kEigenvalueTolerance = 1e-4;
  static const double kEigenvectorTolerance = 1e-4;

  Eigen::MatrixXd matrix(kNumDimensions, kNumDimensions);
  matrix.setRandom();
  matrix = (matrix.transpose() * matrix).eval();
  DenseInverseLULinearOperator linear_operator(matrix);
  DominantEigensolver::Options options;
  DominantEigensolver eigensolver(options, linear_operator);

  // Compute with power iterations.
  Eigen::VectorXd eigenvector;
  double eigenvalue;
  EXPECT_TRUE(eigensolver.Compute(&eigenvalue, &eigenvector));
  // Apply the inverse transformation to recover the proper
  // eigenvalue/eigenvector.
  static const double kTinyNumber = 1e-8;
  eigenvalue = kTinyNumber + 1.0 / eigenvalue;

  // Get ground truth.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gt_eigensolver(matrix);
  const Eigen::VectorXd gt_eigenvector =
      gt_eigensolver.eigenvectors().col(0);
  const double gt_eigenvalue = gt_eigensolver.eigenvalues()(0);

  EXPECT_NEAR(eigenvalue, gt_eigenvalue, kEigenvalueTolerance);
  EXPECT_LT(1.0 - std::abs(eigenvector.dot(gt_eigenvector)),
            kEigenvectorTolerance);
}

TEST(DominantEigensolver, SparseMatrixLargestEigenvalue) {
  static const int kNumDimensions = 100;
  static const double kEigenvalueTolerance = 1e-4;
  static const double kEigenvectorTolerance = 1e-4;

  Eigen::MatrixXd matrix(kNumDimensions, kNumDimensions);
  matrix.setRandom();
  matrix = (matrix.transpose() * matrix).eval();
  const Eigen::SparseMatrix<double> sparse_matrix = matrix.sparseView();
  SparseLinearOperator linear_operator(sparse_matrix);
  DominantEigensolver::Options options;
  options.max_num_iterations = 1000;
  DominantEigensolver eigensolver(options, linear_operator);

  // Compute with power iterations.
  Eigen::VectorXd eigenvector;
  double eigenvalue;
  EXPECT_TRUE(eigensolver.Compute(&eigenvalue, &eigenvector));

  // Get ground truth.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gt_eigensolver(matrix);
  const Eigen::VectorXd gt_eigenvector =
      gt_eigensolver.eigenvectors().rightCols<1>();
  const double gt_eigenvalue = gt_eigensolver.eigenvalues().reverse()(0);

  EXPECT_NEAR(eigenvalue, gt_eigenvalue, kEigenvalueTolerance);
  EXPECT_LT(1.0 - std::abs(eigenvector.dot(gt_eigenvector)),
            kEigenvectorTolerance);
}

TEST(DominantEigensolver, SparseMatrixSmallestEigenvalue) {
  static const int kNumDimensions = 200;
  static const double kEigenvalueTolerance = 1e-4;
  static const double kEigenvectorTolerance = 1e-4;

  Eigen::MatrixXd matrix(kNumDimensions, kNumDimensions);
  matrix.setRandom();
  matrix = (matrix.transpose() * matrix).eval();
  const Eigen::SparseMatrix<double> sparse_matrix = matrix.sparseView();
  SparseInverseLULinearOperator linear_operator(sparse_matrix);
  DominantEigensolver::Options options;
  options.max_num_iterations = 1000;
  DominantEigensolver eigensolver(options, linear_operator);

  // Compute with power iterations.
  Eigen::VectorXd eigenvector;
  double eigenvalue;
  EXPECT_TRUE(eigensolver.Compute(&eigenvalue, &eigenvector));
  // Apply the inverse transformation to recover the proper
  // eigenvalue/eigenvector.
  static const double kTinyNumber = 1e-8;
  eigenvalue = kTinyNumber + 1.0 / eigenvalue;

  // Get ground truth.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gt_eigensolver(matrix);
  const Eigen::VectorXd gt_eigenvector =
      gt_eigensolver.eigenvectors().col(0);
  const double gt_eigenvalue = gt_eigensolver.eigenvalues()(0);
  EXPECT_NEAR(eigenvalue, gt_eigenvalue, kEigenvalueTolerance);
  EXPECT_LT(1.0 - std::abs(eigenvector.dot(gt_eigenvector)),
            kEigenvectorTolerance);
}

TEST(DominantEigensolver, RankDeficiency) {
  static const int kNumDimensions = 500;
  static const double kEigenvalueTolerance = 1e-4;
  static const double kEigenvectorTolerance = 1e-4;

  Eigen::MatrixXd matrix(kNumDimensions, kNumDimensions);
  matrix.setRandom();
  // Force rank deficiency.
  Eigen::VectorXd diag(kNumDimensions);
  diag.setRandom();
  diag = diag.cwiseAbs2();
  diag(0) = 0.0;
  matrix = (matrix.transpose() * diag.asDiagonal() * matrix).eval();

  const Eigen::SparseMatrix<double> sparse_matrix = matrix.sparseView();
  SparseInverseLULinearOperator linear_operator(sparse_matrix);
  DominantEigensolver::Options options;
  DominantEigensolver eigensolver(options, linear_operator);

  // Compute with power iterations.
  Eigen::VectorXd eigenvector;
  double eigenvalue;
  EXPECT_TRUE(eigensolver.Compute(&eigenvalue, &eigenvector));
  // Apply the inverse transformation to recover the proper
  // eigenvalue/eigenvector.
  static const double kTinyNumber = 1e-8;
  eigenvalue = kTinyNumber + 1.0 / eigenvalue;

  // Get ground truth.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gt_eigensolver(matrix);
  const Eigen::VectorXd gt_eigenvector =
      gt_eigensolver.eigenvectors().col(0);
  const double gt_eigenvalue = gt_eigensolver.eigenvalues()(0);

  EXPECT_NEAR(eigenvalue, 0.0, kEigenvalueTolerance);
  EXPECT_NEAR(gt_eigenvalue, 0.0, kEigenvalueTolerance);
  EXPECT_LT(1.0 - std::abs(eigenvector.dot(gt_eigenvector)),
            kEigenvectorTolerance);
}

}  // namespace theia
