// Copyright (C) 2018 The Regents of the University of California (Regents).
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
// Author: Victor Fragoso (victor.fragoso@mail.wvu.edu)

#include "theia/sfm/pose/gauss_jordan_elimination.h"

#include <cmath>
#include <Eigen/Core>
#include <glog/logging.h>

namespace theia {
namespace {
inline double FindLargestValueInColumn(const RowMajorMatrixXd& matrix,
                                       const int column,
                                       const int starting_row,
                                       int* max_value_index) {
  double max_value = matrix(starting_row, column);
  *max_value_index = starting_row;
  for (int row = starting_row; row < matrix.rows(); ++row) {
    if (std::abs(max_value) < std::abs(matrix(row, column))) {
      max_value = matrix(row, column);
      *max_value_index = row;
    }
  }
  return max_value;
}

}  // namespace

// The implementation is based on Gauss Jordan Elimination algorithm explained
// in
// https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
//
// for (row = 0; row < max_num_iterations; ++row) {
//   // Search for the maximum entry in column = row. Keep its value and row.
//   // Swap row with maximum entry to be the current one.
//   // Make all the entries in column zero.
// }
//
// and also inspired by the implementation from OpenGV:
// src/math/gauss_jordan.cpp.
void GaussJordanElimination(const int last_row_to_process,
                            RowMajorMatrixXd* template_matrix) {
  const double kPrecisionThreshold = 1e-9;
  CHECK_GE(last_row_to_process, 0);
  CHECK_LT(last_row_to_process, CHECK_NOTNULL(template_matrix)->rows());
  // Check that the template matrix is a fat matrix.
  CHECK_GE(template_matrix->cols(), template_matrix->rows());

  const int num_rows = template_matrix->rows();
  // This for loop eliminates entries in the lower-left triangular part, and it
  // operates from top to bottom of the matrix.
  for (int current_row = 0; current_row < num_rows; ++current_row) {
    // Search for the maxium entry in column with current_row index.
    const int tail_size = num_rows - current_row;
    int max_coeff_row_idx = 0;
    const double max_value = FindLargestValueInColumn(
        *template_matrix, current_row, current_row, &max_coeff_row_idx);

    // Swap rows.
    template_matrix->row(current_row).
        swap(template_matrix->row(max_coeff_row_idx));

    // Divide the current row by the leading coefficient or max value.
    template_matrix->row(current_row) /= max_value;

    // Set the remaining entries of the column to zero.
    for (int temp_row = current_row + 1; temp_row < num_rows; ++temp_row) {
      // If the number is too small, it is best to consider it as zero already.
      const double leading_coeff = (*template_matrix)(temp_row, current_row);
      if (std::abs(leading_coeff) < kPrecisionThreshold) {
        continue;
      }

      // Multiply current_row by the leading coeff and subtract it from temp_row
      // to eliminate or set the leading coefficient to zero.
      template_matrix->row(temp_row) -=
          leading_coeff * template_matrix->row(current_row);
    }
  }

  // This for loop eliminates entries in the upper-left triangular part, and it
  // operates from bottom to top of the matrix.
  for (int current_row = num_rows - 1;
       current_row >= last_row_to_process; --current_row) {
    // Column to process.
    const int col = current_row;
    // Eliminate coefficients in column.
    for (int row = current_row - 1; row >= last_row_to_process; --row) {
      // Get the leading coefficient to eliminate.
      const double leading_coeff = (*template_matrix)(row, col);
      if (std::abs(leading_coeff) < kPrecisionThreshold) {
        continue;
      }

      // Eliminate.
      template_matrix->row(row) -=
          leading_coeff * template_matrix->row(current_row);
    }
  }
}

}  // namespace theia
