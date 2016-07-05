//=============================================================================
//
// nldiffusion_functions.cpp
// Authors: Pablo F. Alcantarilla (1), Jesus Nuevo (2)
// Institutions: Toshiba Research Europe Ltd (1)
//               TrueVision Solutions (2)
// Date: 07/10/2014
// Email: pablofdezalc@gmail.com
//
// AKAZE Features Copyright 2014, Pablo F. Alcantarilla, Jesus Nuevo
// All Rights Reserved
// See LICENSE for the license information
//=============================================================================

/**
 * @file nldiffusion_functions.cpp
 * @brief Functions for nonlinear diffusion filtering applications
 * @date Oct 07, 2014
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#include <cassert>
#ifdef AKAZE_USE_OPENMP
#include <omp.h>
#endif  // AKAZE_USE_OPENMP
#include <vector>

#include "convolution.h"
#include "nldiffusion_functions.h"

namespace libAKAZE {

/* ************************************************************************* */
void gaussian_2D_convolution(const RowMatrixXf& src, RowMatrixXf& dst,
                             size_t ksize_x, size_t ksize_y, float sigma) {
  GaussianBlur(src, sigma, &dst);
}

/* ************************************************************************* */
void image_derivatives_scharr(const RowMatrixXf& src, RowMatrixXf& dst,
                              const size_t xorder, const size_t yorder) {
  ScharrDerivative(src, xorder, yorder, 1.0, false, &dst);
}

/* ************************************************************************* */
void pm_g1(const RowMatrixXf& Lx, const RowMatrixXf& Ly, RowMatrixXf& dst,
           const float k) {
  const float inv_k = 1.0 / (k * k);
  dst = (-inv_k * (Lx.array().square() + Ly.array().square())).exp();
}

/* ************************************************************************* */
void pm_g2(const RowMatrixXf& Lx, const RowMatrixXf& Ly, RowMatrixXf& dst,
           const float k) {
  const float inv_k = 1.0 / (k * k);
  dst = (1.0 + inv_k * (Lx.array().square() + Ly.array().square())).inverse();
}

/* ************************************************************************* */
void weickert_diffusivity(const RowMatrixXf& Lx, const RowMatrixXf& Ly,
                          RowMatrixXf& dst, const float k) {
  dst.resize(Lx.rows(), Lx.cols());
  const float inv_k = 1.0 / (k * k);
  for (int y = 0; y < Lx.rows(); y++) {
    for (int x = 0; x < Lx.cols(); x++) {
      const float dL = inv_k * (Lx(y, x) * Lx(y, x) + Ly(y, x) * Ly(y, x));
      dst(y, x) = -3.315 / (dL * dL * dL * dL);
    }
  }
  dst = dst.array().exp();
  dst = 1.0 - dst.array();
}

/* ************************************************************************* */
void charbonnier_diffusivity(const RowMatrixXf& Lx, const RowMatrixXf& Ly,
                             RowMatrixXf& dst, const float k) {
  const float inv_k = 1.0 / (k * k);
  dst = 1.0 /
        ((1.0 + inv_k * (Lx.array().square() + Ly.array().square())).sqrt());
}

/* ************************************************************************* */
float compute_k_percentile(const RowMatrixXf& img, float perc, float gscale,
                           size_t nbins, size_t ksize_x, size_t ksize_y) {
  size_t nbin = 0, nelements = 0, nthreshold = 0, k = 0;
  float kperc = 0.0, modg = 0.0, npoints = 0.0, hmax = 0.0;

  // Create the array for the histogram
  Eigen::VectorXf hist(nbins);
  hist.setZero();

  // Create the matrices
  RowMatrixXf gaussian(img.rows(), img.cols());
  RowMatrixXf Lx(img.rows(), img.cols());
  RowMatrixXf Ly(img.rows(), img.cols());

  // Perform the Gaussian convolution
  gaussian_2D_convolution(img, gaussian, ksize_x, ksize_y, gscale);

  // Compute the Gaussian derivatives Lx and Ly
  image_derivatives_scharr(gaussian, Lx, 1, 0);
  image_derivatives_scharr(gaussian, Ly, 0, 1);

  // Skip the borders for computing the histogram
  for (int y = 1; y < gaussian.rows() - 1; y++) {
    for (int x = 1; x < gaussian.cols() - 1; x++) {
      modg = Lx(y, x) * Lx(y, x) + Ly(y, x) * Ly(y, x);
      // Get the maximum
      if (modg > hmax) hmax = modg;
    }
  }

  hmax = sqrt(hmax);

  // Skip the borders for computing the histogram
  for (int y = 1; y < gaussian.rows() - 1; y++) {
    for (int x = 1; x < gaussian.cols() - 1; x++) {
      modg = sqrt(Lx(y, x) * Lx(y, x) + Ly(y, x) * Ly(y, x));

      // Find the correspondent bin
      if (modg != 0.0) {
        nbin = floor(nbins * (modg / hmax));

        if (nbin == nbins) {
          nbin--;
        }

        hist[nbin]++;
        npoints++;
      }
    }
  }

  // Now find the perc of the histogram percentile
  nthreshold = (size_t)(npoints * perc);

  for (k = 0; nelements < nthreshold && k < nbins; k++)
    nelements = nelements + hist[k];

  if (nelements < nthreshold)
    kperc = 0.03;
  else
    kperc = hmax * ((float)(k) / (float) nbins);

  return kperc;
}

/* ************************************************************************* */
void compute_scharr_derivatives(const RowMatrixXf& src,
                                RowMatrixXf& dst,
                                const size_t xorder, const size_t yorder,
                                const size_t scale) {
  ScharrDerivative(src, xorder, yorder, scale, true, &dst);
}

/* ************************************************************************* */
void nld_step_scalar(RowMatrixXf& Ld, const RowMatrixXf& c, RowMatrixXf& Lstep,
                     const float stepsize) {
  Lstep.resize(Ld.rows(), Ld.cols());
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int y = 1; y < Lstep.rows() - 1; y++) {
    for (int x = 1; x < Lstep.cols() - 1; x++) {
      const float xpos = (c(y, x) + c(y, x + 1)) * (Ld(y, x + 1) - Ld(y, x));
      const float xneg = (c(y, x - 1) + c(y, x)) * (Ld(y, x) - Ld(y, x - 1));
      const float ypos = (c(y, x) + c(y + 1, x)) * (Ld(y + 1, x) - Ld(y, x));
      const float yneg = (c(y - 1, x) + c(y, x)) * (Ld(y, x) - Ld(y - 1, x));
      Lstep(y, x) = 0.5 * stepsize * (xpos - xneg + ypos - yneg);
    }
  }

  for (int x = 1; x < Lstep.cols() - 1; x++) {
    const float xpos = (c(0, x) + c(0, x + 1)) * (Ld(0, x + 1) - Ld(0, x));
    const float xneg = (c(0, x - 1) + c(0, x)) * (Ld(0, x) - Ld(0, x - 1));
    const float ypos = (c(0, x) + c(1, x)) * (Ld(1, x) - Ld(0, x));
    Lstep(0, x) = 0.5 * stepsize * (xpos - xneg + ypos);
  }

  const int end_row = Lstep.rows() - 1;
  for (int x = 1; x < Lstep.cols() - 1; x++) {
    const float xpos = (c(end_row, x) + c(end_row, x + 1)) *
                       (Ld(end_row, x + 1) - Ld(end_row, x));
    const float xneg = (c(end_row, x - 1) + c(end_row, x)) *
                       (Ld(end_row, x) - Ld(end_row, x - 1));
    const float ypos = (c(end_row, x) + c(end_row - 1, x)) *
                       (Ld(end_row - 1, x) - Ld(end_row, x));
    Lstep(end_row, x) = 0.5 * stepsize * (xpos - xneg + ypos);
  }

  const int c_min_1 = Lstep.cols() - 1;
  const int c_min_2 = Lstep.cols() - 2;
  for (int i = 1; i < Lstep.rows() - 1; i++) {
    float xpos = (c(i, 0) + c(i, 1)) * (Ld(i, 1) - Ld(i, 0));
    float ypos = (c(i, 0) + c(i + 1, 0)) * (Ld(i + 1, 0) - Ld(i, 0));
    float yneg = (c(i - 1, 0) + c(i, 0)) * (Ld(i, 0) - Ld(i - 1, 0));
    Lstep(i, 0) = 0.5 * stepsize * (xpos + ypos - yneg);

    const float xneg =
        (c(i, c_min_2) + c(i, c_min_1)) * (Ld(i, c_min_1) - Ld(i, c_min_2));
    ypos = (c(i, c_min_1) + c(i + 1, c_min_1)) *
           (Ld(i + 1, c_min_1) - Ld(i, c_min_1));
    yneg = (c(i - 1, c_min_1) + c(i, c_min_1)) *
           (Ld(i, c_min_1) - Ld(i - 1, c_min_1));
    Lstep(i, c_min_1) = 0.5 * stepsize * (-xneg + ypos - yneg);
  }

  // Ld = Ld + Lstep
  Ld += Lstep;
}

/* ************************************************************************* */
// I think OpenCV's inter area works by setting the size of the interpolation
// kernel to be equal to the scaling factor.
//
// TODO: OpenCV is ~7x faster for this method when dimensions are odd. Maybe the
// difference is only in using OpenMP (I tested it on a Mac). Should try to
// improve the performance here.
void halfsample_image(const RowMatrixXf& src, RowMatrixXf& dst) {
  assert(src.rows() / 2 == dst.rows());
  assert(src.cols() / 2 == dst.cols());

  // Do the whole resampling in one pass by using neighboring values. First, we
  // compute the borders.
  const double x_kernel_size = static_cast<double>(src.cols()) / dst.cols();
  const double y_kernel_size = static_cast<double>(src.rows()) / dst.rows();

  // Do simple linear interpolation.
  if (x_kernel_size == 2 && y_kernel_size == 2) {
    for (int i = 0; i < dst.rows(); i++) {
      for (int j = 0; j < dst.cols(); j++) {
        dst(i, j) = src(2 * i, 2 * j) + src(2 * i + 1, 2 * j) +
                    src(2 * i, 2 * j + 1) + src(2 * i + 1, 2 * j + 1);
      }
    }
    dst /= 4.0;
    return;
  }

  const double x_kernel_clamped_size = static_cast<int>(ceil(x_kernel_size));
  const double y_kernel_clamped_size = static_cast<int>(ceil(y_kernel_size));

  // Set up precomputed factor matrices.
  Eigen::RowVectorXf x_kernel_mul(static_cast<int>(x_kernel_clamped_size)),
      y_kernel_mul(static_cast<int>(y_kernel_clamped_size));
  y_kernel_mul.setConstant(1.0);

  Eigen::RowVectorXf temp_row(src.cols());
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for firstprivate(temp_row, x_kernel_mul, y_kernel_mul)
#endif
  for (int i = 0; i < dst.rows(); i++) {
    // Compute the row resize first.
    const int y = static_cast<int>(y_kernel_size * i);
    y_kernel_mul(0) = 1 - (y_kernel_size * i - y);
    y_kernel_mul(y_kernel_clamped_size - 1) -=
        y_kernel_mul.sum() - y_kernel_size;

    temp_row =
        y_kernel_mul * src.block(y, 0, y_kernel_clamped_size, src.cols());

    // For this row, compute the column-wise resize.
    x_kernel_mul.setConstant(1.0);
    for (int j = 0; j < dst.cols(); j++) {
      const int x = static_cast<int>(x_kernel_size * j);
      x_kernel_mul(0) = 1 - (x_kernel_size * j - x);
      x_kernel_mul(x_kernel_clamped_size - 1) -=
          x_kernel_mul.sum() - x_kernel_size;
      dst(i, j) =
          x_kernel_mul.dot(temp_row.segment(x, x_kernel_clamped_size));
    }
  }

  dst /= x_kernel_size * y_kernel_size;
}

/* ************************************************************************* */
bool check_maximum_neighbourhood(const RowMatrixXf& img, int dsize, float value,
                                 int row, int col, bool same_img) {
  bool response = true;

  for (int i = row - dsize; i <= row + dsize; i++) {
    for (int j = col - dsize; j <= col + dsize; j++) {
      if (i >= 0 && i < img.rows() && j >= 0 && j < img.cols()) {
        if (same_img == true) {
          if (i != row || j != col) {
            if (img(i, j) > value) {
              response = false;
              return response;
            }
          }
        } else {
          if (img(i, j) > value) {
            response = false;
            return response;
          }
        }
      }
    }
  }

  return response;
}

}  // namespace libAKAZE
