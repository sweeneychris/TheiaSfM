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

#ifndef THEIA_MATH_MATRIX_LINEAR_OPERATOR_H_
#define THEIA_MATH_MATRIX_LINEAR_OPERATOR_H_

#include <Eigen/SparseCore>
#include <glog/logging.h>

#include "theia/math/matrix/sparse_cholesky_llt.h"

namespace theia {

// A sparse method for computing the shift and inverse linear operator. This
// method is intended for use with the Spectra library.
class SparseSymShiftSolveLLT {
 public:
  explicit SparseSymShiftSolveLLT(const Eigen::SparseMatrix<double>& mat)
      : mat_(mat) {
    CHECK_EQ(mat_.rows(), mat_.cols());
    linear_solver_.Compute(mat_);
    if (linear_solver_.Info() != Eigen::Success) {
      LOG(FATAL)
          << "Could not perform Cholesky decomposition on the matrix. Are "
          "you sure it is positive semi-definite?";
    }
  }

  int rows() { return mat_.rows(); }
  int cols() { return mat_.cols(); }
  void set_shift(double sigma) { sigma_ = sigma; }

  // Use LDLT to perform matrix inversion on the positive semidefinite matrix.
  void perform_op(double* x_in, double* y_out) {
    Eigen::Map<Eigen::VectorXd> x(x_in, mat_.rows());
    Eigen::Map<Eigen::VectorXd> y(y_out, mat_.cols());
    y = linear_solver_.Solve(x);
    if (linear_solver_.Info() != Eigen::Success) {
      LOG(FATAL)
          << "Could not perform Cholesky decomposition on the matrix. Are "
          "you sure it is positive semi-definite?";
    }
  }

  const Eigen::SparseMatrix<double>& mat_;
  SparseCholeskyLLt linear_solver_;
  double sigma_;
};

}  // namespace theia

#endif  // THEIA_MATH_MATRIX_LINEAR_OPERATOR_H_
