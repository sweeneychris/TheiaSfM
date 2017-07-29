// Copyright (C) 2017 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#ifndef THEIA_MATH_MATRIX_SPARSE_MATRIX_H_
#define THEIA_MATH_MATRIX_SPARSE_MATRIX_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "theia/math/matrix/linear_operator.h"

namespace theia {

// A wrapper for the Eigen sparse matrix class.
class SparseMatrix : public LinearOperator {
 public:
  explicit SparseMatrix(const Eigen::SparseMatrix<double>& matrix);
  ~SparseMatrix();

  // y = y + Ax;
  void RightMultiply(const Eigen::VectorXd& x,
                     Eigen::VectorXd* y) const override;

  // y = y + A'x;
  void LeftMultiply(const Eigen::VectorXd& x,
                    Eigen::VectorXd* y) const override;

  int num_rows() const override;
  int num_cols() const override;

 private:
  const Eigen::SparseMatrix<double>& matrix_;
};  // namespace theia

}  // namespace theia

#endif  // THEIA_MATH_MATRIX_SPARSE_MATRIX_H_
