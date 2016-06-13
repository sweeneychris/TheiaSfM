// Copyright (C) 2014 The Regents of the University of California (Regents).
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

#include "convolution.h"

#include <Eigen/Dense>
#include <iostream>

#ifdef AKAZE_USE_OPENMP
#include <omp.h>
#endif  // AKAZE_USE_OPENMP

namespace {

using Eigen::Matrix;
using Eigen::RowVectorXf;

inline double Gaussian(double x, double mu, double sigma) {
  return std::exp(-(x - mu) * (x - mu) / (2.0 * sigma * sigma));
}

}  // namespace

void SeparableConvolution2d(const RowMatrixXf& image,
                            const Eigen::RowVectorXf& kernel_x,
                            const Eigen::RowVectorXf& kernel_y,
                            const BorderType& border_type,
                            RowMatrixXf* out) {
  const int full_size = kernel_x.size();
  const int half_size = full_size / 2;
  out->resize(image.rows(), image.cols());

  // Convolving a vertical filter across rows is the same thing as transpose
  // multiply i.e. kernel_y^t * rows. This will give us the convoled value for
  // each row. However, care must be taken at the top and bottom borders.
  const RowVectorXf reverse_kernel_y = kernel_y.reverse();

  if (border_type == REFLECT) {
    for (int i = 0; i < half_size; i++) {
      const int forward_size = i + half_size + 1;
      const int reverse_size = full_size - forward_size;
      out->row(i) = kernel_y.tail(forward_size) *
                    image.block(0, 0, forward_size, image.cols()) +
                    reverse_kernel_y.tail(reverse_size) *
                    image.block(1, 0, reverse_size, image.cols());

      // Apply the same technique for the end rows.
      // TODO(csweeney): Move this to its own loop for cache exposure?
      out->row(image.rows() - i - 1) =
          kernel_y.head(forward_size) * image.block(image.rows() - forward_size,
                                                    0, forward_size,
                                                    image.cols()) +
          reverse_kernel_y.head(reverse_size) *
          image.block(image.rows() - reverse_size - 1, 0, reverse_size,
                      image.cols());
    }
  } else {
    // Perform border with REPLICATE as the option.
    for (int i = 0; i < half_size; i++) {
      const int forward_size = i + half_size + 1;
      const int reverse_size = full_size - forward_size;
      out->row(i) = kernel_y.tail(forward_size) *
                        image.block(0, 0, forward_size, image.cols()) +
                    reverse_kernel_y.tail(reverse_size) *
                        image.row(0).replicate(reverse_size, 1);

      // Apply the same technique for the end rows.
      out->row(image.rows() - i - 1) =
          kernel_y.head(forward_size) * image.block(image.rows() - forward_size,
                                                    0, forward_size,
                                                    image.cols()) +
          reverse_kernel_y.head(reverse_size) *
          image.row(image.rows() - 1).replicate(reverse_size, 1);
    }
  }

  // Applying the rest of the y filter.
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int row = half_size; row < image.rows() - half_size; row++) {
    out->row(row) =
        kernel_y * image.block(row - half_size, 0, full_size, out->cols());
  }

  // Convolving with the horizontal filter is easy. Rather than using the kernel
  // as a sliding indow, we use the row pixels as a sliding window around the
  // filter. We prepend and append the proper border values so that we are sure
  // to end up with the correct convolved values.
  if (border_type == REFLECT) {
    RowVectorXf temp_row(image.cols() + full_size - 1);
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for firstprivate(temp_row)
#endif
    for (int row = 0; row < out->rows(); row++) {
    temp_row.head(half_size) =
          out->row(row).segment(1, half_size).reverse();
      temp_row.segment(half_size, image.cols()) = out->row(row);
      temp_row.tail(half_size) =
          out->row(row)
          .segment(image.cols() - 1 - half_size, half_size)
          .reverse();

      // Convolve the row. We perform the first step here explicitly so that we
      // avoid setting the row equal to zero.
      out->row(row) = kernel_x(0) * temp_row.head(image.cols());
      for (int i = 1; i < full_size; i++) {
        out->row(row) += kernel_x(i) * temp_row.segment(i, image.cols());
      }
    }
  } else {
    RowVectorXf temp_row(image.cols() + full_size - 1);
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for firstprivate(temp_row)
#endif
    for (int row = 0; row < out->rows(); row++) {
      temp_row.head(half_size).setConstant((*out)(row, 0));
      temp_row.segment(half_size, image.cols()) = out->row(row);
      temp_row.tail(half_size).setConstant((*out)(row, out->cols() - 1));

      // Convolve the row. We perform the first step here explicitly so that we
      // avoid setting the row equal to zero.
      out->row(row) = kernel_x(0) * temp_row.head(image.cols());
      for (int i = 1; i < full_size; i++) {
        out->row(row) += kernel_x(i) * temp_row.segment(i, image.cols());
      }
    }
  }
}

void ScharrDerivative(const RowMatrixXf& image,
                      const int x_deg,
                      const int y_deg,
                      const int size,
                      const bool normalize,
                      RowMatrixXf* out) {
  const int sigma = size * 2 + 1;
  Eigen::RowVectorXf kernel1(sigma);
  kernel1.setZero();
  kernel1(0) = -1.0;
  kernel1(sigma - 1) = 1.0;

  Eigen::RowVectorXf kernel2(sigma);
  kernel2.setZero();
  if (!normalize) {
    kernel2(0) = 3;
    kernel2(sigma / 2) = 10;
    kernel2(sigma - 1) = 3;
  } else {
    float w = 10.0 / 3.0;
    float norm = 1.0 / (2.0 * size * (w + 2.0));
    kernel2(0) = norm;
    kernel2(sigma / 2) = w * norm;
    kernel2(sigma - 1) = norm;
  }

  if (x_deg == 1) {
    SeparableConvolution2d(image, kernel1, kernel2, REFLECT, out);
  } else {
    SeparableConvolution2d(image, kernel2, kernel1, REFLECT, out);
  }
  return;
}

void GaussianBlur(const RowMatrixXf& image,
                  const double sigma,
                  RowMatrixXf* out) {
  int kernel_size = std::ceil(((sigma - 0.8) / 0.3 + 1.0) * 2.0);
  if (kernel_size % 2 == 0) {
    kernel_size += 1;
  }

  RowVectorXf gauss_kernel(kernel_size);
  double norm_factor = 0;
  for (int i = 0; i < gauss_kernel.size(); i++) {
    gauss_kernel(i) = Gaussian(i, (kernel_size - 1.0) / 2.0, sigma);
    norm_factor += gauss_kernel(i);
  }
  gauss_kernel /= norm_factor;

  SeparableConvolution2d(image,
                         gauss_kernel,
                         gauss_kernel,
                         REPLICATE,
                         out);
}
