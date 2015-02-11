// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#ifndef THEIA_MATH_MATRIX_GAUSS_JORDAN_H_
#define THEIA_MATH_MATRIX_GAUSS_JORDAN_H_

#include <Eigen/Core>
#include <algorithm>

namespace theia {

// Gauss-Jordan elimination on a matrix.
// Modifies the input matrix to be the matrix after gauss-jordan elimation.
template <typename Derived>
void GaussJordan(Eigen::MatrixBase<Derived>* input, int max_rows = 99999) {
  max_rows = std::min(static_cast<int>(input->rows()), max_rows);
  // Iterate through each column one by one.
  for (int i = 0; i < max_rows; i++) {
    // Find row with the largest value in the pivot column and swap.
    int swap_row = i;
    double max_val = 0.0;
    for (int j = i + 1; j < input->rows(); j++) {
      double temp_max_val = std::abs((*input)(j, i));
      if (temp_max_val > max_val) {
        max_val = temp_max_val;
        swap_row = j;
      }
    }
    input->row(swap_row).swap(input->row(i));

    // Eliminate all values in the column of the pivot.
    for (int j = 0; j < input->rows(); j++) {
      if (j != i) {
        input->row(j) -= ((*input)(j, i) / (*input)(i, i)) * input->row(i);
        (*input)(j, i) = 0.0;
      }
    }
  }

  for (int i = 0; i < max_rows; i++) {
    // Scale current row so that leading value is 1.0. Leading value should be
    // along the diagonal as we proceed with gauss-jordan.
    input->row(i) *= 1.0 / (*input)(i, i);
    (*input)(i, i) = 1.0;
  }
}
}  // namespace theia

#endif  // THEIA_MATH_MATRIX_GAUSS_JORDAN_H_
