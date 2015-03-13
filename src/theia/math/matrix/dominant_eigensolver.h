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

#ifndef THEIA_MATH_MATRIX_DOMINANT_EIGENSOLVER_H_
#define THEIA_MATH_MATRIX_DOMINANT_EIGENSOLVER_H_

#include <Eigen/Core>

#include "theia/math/matrix/linear_operator.h"

namespace theia {

// Computes the dominant eigenvalue of the matrix specified by the linear
// operator A by using the Power method. This method can also be used to compute
// the smallest eigenvalue or the eigenvalue nearest to a value b by
// implementing (A - b * I)^-1 as the RightMultiply of the linear operator. This
// can be done efficiently with a linear solve.
class DominantEigensolver {
 public:
  struct Options {
    // Maximum number of iterations to run.
    int max_num_iterations = 100;

    // The tolerance for determining when a solution has converged.
    double tolerance = 1e-6;
  };

  DominantEigensolver(const Options& options, const LinearOperator& A)
      : options_(options), A_(A) {}

  // Computes the dominant eigenvalue and eigenvector using the Power iteration
  // method.
  bool Compute(double* eigenvalue, Eigen::VectorXd* eigenvector) const;

 private:
  const Options options_;
  const LinearOperator& A_;
};

}  // namespace theia

#endif  // THEIA_MATH_MATRIX_DOMINANT_EIGENSOLVER_H_
