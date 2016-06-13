//=============================================================================
//
// AKAZE.cpp
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
 * @file AKAZE.cpp
 * @brief Main class for detecting and describing binary features in an
 * accelerated nonlinear scale space
 * @date Oct 07, 2014
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#include "cimg/CImg.h"
#include "src/AKAZE.h"
#include "src/utils.h"
#include "timer/timer.hpp"

#include <Eigen/Dense>
#include <math.h>
#include <cassert>
#include <ctime>
#include <iostream>

#define _USE_MATH_DEFINESo

using namespace libAKAZE;

/* ************************************************************************* */
AKAZE::AKAZE(const AKAZEOptions& options) : options_(options) {
  Eigen::initParallel();
#ifdef AKAZE_USE_OPENMP
  omp_set_num_threads(options.num_threads);
#define EIGEN_DONT_PARALLELIZE
#endif
  ncycles_ = 0;
  reordering_ = true;

  if (options_.descriptor_size > 0 && options_.descriptor >= MLDB_UPRIGHT) {
    generateDescriptorSubsample(
        descriptorSamples_, descriptorBits_, options_.descriptor_size,
        options_.descriptor_pattern_size, options_.descriptor_channels);
  }
  Allocate_Memory_Evolution();
}

/* ************************************************************************* */
AKAZE::~AKAZE() { evolution_.clear(); }

/* ************************************************************************* */
void AKAZE::Allocate_Memory_Evolution() {

  float rfactor = 0.0;
  int level_height = 0, level_width = 0;

  // Allocate the dimension of the matrices for the evolution
  evolution_.reserve(options_.omax * options_.nsublevels);
  for (int i = 0; i < options_.omax; i++) {
    rfactor = 1.0 / pow(2.f, i);
    level_height = (int)(options_.img_height * rfactor);
    level_width = (int)(options_.img_width * rfactor);

    // Smallest possible octave and allow one scale if the image is small
    if ((level_width < 80 || level_height < 40) && i != 0) {
      options_.omax = i;
      break;
    }

    for (int j = 0; j < options_.nsublevels; j++) {
      TEvolution step;
      step.Lx.resize(level_height, level_width);
      step.Ly.resize(level_height, level_width);
      step.Lxx.resize(level_height, level_width);
      step.Lxy.resize(level_height, level_width);
      step.Lyy.resize(level_height, level_width);
      step.Lt.resize(level_height, level_width);
      step.Ldet.resize(level_height, level_width);
      step.Lflow.resize(level_height, level_width);
      step.Lstep.resize(level_height, level_width);

      step.esigma = options_.soffset *
                    pow(2.f, (float)(j) / (float)(options_.nsublevels) + i);
      step.sigma_size = fRound(step.esigma);
      step.etime = 0.5 * (step.esigma * step.esigma);
      step.octave = i;
      step.sublevel = j;
      evolution_.push_back(step);
    }
  }

  // Allocate memory for the number of cycles and time steps
  for (size_t i = 1; i < evolution_.size(); i++) {
    int naux = 0;
    std::vector<float> tau;
    float ttime = 0.0;
    ttime = evolution_[i].etime - evolution_[i - 1].etime;
    naux = fed_tau_by_process_time(ttime, 1, 0.25, reordering_, tau);
    nsteps_.push_back(naux);
    tsteps_.push_back(tau);
    ncycles_++;
  }
}

/* ************************************************************************* */
int AKAZE::Create_Nonlinear_Scale_Space(const RowMatrixXf& img) {
  if (evolution_.size() == 0) {
    std::cerr << "Error generating the nonlinear scale space!!" << std::endl;
    std::cerr << "Firstly you need to call AKAZE::Allocate_Memory_Evolution()"
              << std::endl;
    return -1;
  }

  timer::Timer timer;

  // Copy the original image to the first level of the evolution
  gaussian_2D_convolution(img, evolution_[0].Lt, 0, 0, options_.soffset);
  evolution_[0].Lsmooth = evolution_[0].Lt;

  // First compute the kcontrast factor
  options_.kcontrast = compute_k_percentile(
      img, options_.kcontrast_percentile, 1.0, options_.kcontrast_nbins, 0, 0);

  timing_.kcontrast = timer.elapsedMs();

  // Now generate the rest of evolution levels
  for (size_t i = 1; i < evolution_.size(); i++) {

    if (evolution_[i].octave > evolution_[i - 1].octave) {
      halfsample_image(evolution_[i - 1].Lt, evolution_[i].Lt);
      options_.kcontrast = options_.kcontrast * 0.75;
    } else {
      evolution_[i].Lt = evolution_[i - 1].Lt;
    }

    gaussian_2D_convolution(evolution_[i].Lt, evolution_[i].Lsmooth, 0, 0, 1.0);

    // Compute the Gaussian derivatives Lx and Ly
    image_derivatives_scharr(evolution_[i].Lsmooth, evolution_[i].Lx, 1, 0);
    image_derivatives_scharr(evolution_[i].Lsmooth, evolution_[i].Ly, 0, 1);

    // Compute the conductivity equation
    switch (options_.diffusivity) {
      case PM_G1:
        pm_g1(evolution_[i].Lx, evolution_[i].Ly, evolution_[i].Lflow,
              options_.kcontrast);
        break;
      case PM_G2:
        pm_g2(evolution_[i].Lx, evolution_[i].Ly, evolution_[i].Lflow,
              options_.kcontrast);
        break;
      case WEICKERT:
        weickert_diffusivity(evolution_[i].Lx, evolution_[i].Ly,
                             evolution_[i].Lflow, options_.kcontrast);
        break;
      case CHARBONNIER:
        charbonnier_diffusivity(evolution_[i].Lx, evolution_[i].Ly,
                                evolution_[i].Lflow, options_.kcontrast);
        break;
      default:
        std::cerr << "Diffusivity: " << options_.diffusivity
                  << " is not supported" << std::endl;
    }

    // Perform FED n inner steps
    for (int j = 0; j < nsteps_[i - 1]; j++)
      nld_step_scalar(evolution_[i].Lt, evolution_[i].Lflow,
                      evolution_[i].Lstep, tsteps_[i - 1][j]);
  }

  timing_.scale = timer.elapsedMs();

  return 0;
}

/* ************************************************************************* */
void AKAZE::Feature_Detection(std::vector<AKAZEKeypoint>& kpts) {
  timer::Timer timer;

  kpts.clear();
  Compute_Determinant_Hessian_Response();
  Find_Scale_Space_Extrema(kpts);
  Do_Subpixel_Refinement(kpts);

  timing_.detector = timer.elapsedMs();
}

/* ************************************************************************* */
void AKAZE::Compute_Multiscale_Derivatives() {
  timer::Timer timer;

#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < (int)(evolution_.size()); i++) {

    float ratio = pow(2.f, (float)evolution_[i].octave);
    int sigma_size_ =
        fRound(evolution_[i].esigma * options_.derivative_factor / ratio);

    compute_scharr_derivatives(evolution_[i].Lsmooth, evolution_[i].Lx, 1, 0,
                               sigma_size_);
    compute_scharr_derivatives(evolution_[i].Lsmooth, evolution_[i].Ly, 0, 1,
                               sigma_size_);
    compute_scharr_derivatives(evolution_[i].Lx, evolution_[i].Lxx, 1, 0,
                               sigma_size_);
    compute_scharr_derivatives(evolution_[i].Ly, evolution_[i].Lyy, 0, 1,
                               sigma_size_);
    compute_scharr_derivatives(evolution_[i].Lx, evolution_[i].Lxy, 0, 1,
                               sigma_size_);

    evolution_[i].Lx = evolution_[i].Lx * ((sigma_size_));
    evolution_[i].Ly = evolution_[i].Ly * ((sigma_size_));
    evolution_[i].Lxx = evolution_[i].Lxx * ((sigma_size_) * (sigma_size_));
    evolution_[i].Lxy = evolution_[i].Lxy * ((sigma_size_) * (sigma_size_));
    evolution_[i].Lyy = evolution_[i].Lyy * ((sigma_size_) * (sigma_size_));
  }

  timing_.derivatives = timer.elapsedMs();
}

/* ************************************************************************* */
void AKAZE::Compute_Determinant_Hessian_Response() {

  // Firstly compute the multiscale derivatives
  Compute_Multiscale_Derivatives();

  for (size_t i = 0; i < evolution_.size(); i++) {
    if (options_.verbosity == true) {
      std::cout
          << "Computing detector response. Determinant of Hessian. Evolution "
             "time: " << evolution_[i].etime << std::endl;
    }

    evolution_[i].Ldet =
        evolution_[i].Lxx.cwiseProduct(evolution_[i].Lyy).array() -
        evolution_[i].Lxy.array().square();
  }
}

/* ************************************************************************* */
void AKAZE::Find_Scale_Space_Extrema(std::vector<AKAZEKeypoint>& kpts) {
  float value = 0.0;
  float dist = 0.0, ratio = 0.0, smax = 0.0;
  int npoints = 0, id_repeated = 0;
  int sigma_size_ = 0, left_x = 0, right_x = 0, up_y = 0, down_y = 0;
  bool is_extremum = false, is_repeated = false, is_out = false;
  AKAZEKeypoint point;
  std::vector<AKAZEKeypoint> kpts_aux;

  // Set maximum size
  if (options_.descriptor == SURF_UPRIGHT || options_.descriptor == SURF ||
      options_.descriptor == MLDB_UPRIGHT || options_.descriptor == MLDB) {
    smax = 10.0 * sqrtf(2.0);
  } else if (options_.descriptor == MSURF_UPRIGHT ||
             options_.descriptor == MSURF) {
    smax = 12.0 * sqrtf(2.0);
  }

  timer::Timer timer;

  for (size_t i = 0; i < evolution_.size(); i++) {
    for (int ix = 1; ix < evolution_[i].Ldet.rows() - 1; ix++) {
      for (int jx = 1; jx < evolution_[i].Ldet.cols() - 1; jx++) {
        is_extremum = false;
        is_repeated = false;
        is_out = false;
        RowMatrixXf& ldet = evolution_[i].Ldet;
        value = ldet(ix, jx);

        // Filter the points with the detector threshold
        if (value > options_.dthreshold && value >= options_.min_dthreshold &&
            value > ldet(ix, jx - 1) && value > ldet(ix, jx + 1) &&
            value > ldet(ix - 1, jx - 1) && value > ldet(ix - 1, jx) &&
            value > ldet(ix - 1, jx + 1) && value > ldet(ix + 1, jx - 1) &&
            value > ldet(ix + 1, jx) && value > ldet(ix + 1, jx + 1)) {

          is_extremum = true;
          point.response = fabs(value);
          point.size = evolution_[i].esigma * options_.derivative_factor;
          point.octave = evolution_[i].octave;
          point.class_id = i;
          ratio = pow(2.f, point.octave);
          sigma_size_ = fRound(point.size / ratio);
          point.pt.x() = jx;
          point.pt.y() = ix;

          // Compare response with the same and lower scale
          for (size_t ik = 0; ik < kpts_aux.size(); ik++) {

            if ((point.class_id - 1) == kpts_aux[ik].class_id ||
                point.class_id == kpts_aux[ik].class_id) {
              dist = sqrt(pow(point.pt.x() * ratio - kpts_aux[ik].pt.x(), 2) +
                          pow(point.pt.y() * ratio - kpts_aux[ik].pt.y(), 2));
              if (dist <= point.size) {
                if (point.response > kpts_aux[ik].response) {
                  id_repeated = ik;
                  is_repeated = true;
                } else {
                  is_extremum = false;
                }
                break;
              }
            }
          }

          // Check out of bounds
          if (is_extremum == true) {

            // Check that the point is under the image limits for the descriptor
            // computation
            left_x = fRound(point.pt.x() - smax * sigma_size_) - 1;
            right_x = fRound(point.pt.x() + smax * sigma_size_) + 1;
            up_y = fRound(point.pt.y() - smax * sigma_size_) - 1;
            down_y = fRound(point.pt.y() + smax * sigma_size_) + 1;

            if (left_x < 0 || right_x >= evolution_[i].Ldet.cols() ||
                up_y < 0 || down_y >= evolution_[i].Ldet.rows()) {
              is_out = true;
            }

            if (is_out == false) {
              if (is_repeated == false) {
                point.pt.x() *= ratio;
                point.pt.y() *= ratio;
                kpts_aux.push_back(point);
                npoints++;
              } else {
                point.pt.x() *= ratio;
                point.pt.y() *= ratio;
                kpts_aux[id_repeated] = point;
              }
            }  // if is_out
          }    // if is_extremum
        }
      }  // for jx
    }    // for ix
  }      // for i

  // Now filter points with the upper scale level
  kpts.reserve(kpts_aux.size());
  for (size_t i = 0; i < kpts_aux.size(); i++) {

    is_repeated = false;
    const AKAZEKeypoint& point = kpts_aux[i];
    for (size_t j = i + 1; j < kpts_aux.size(); j++) {

      // Compare response with the upper scale
      if ((point.class_id + 1) == kpts_aux[j].class_id) {
        dist = (point.pt - kpts_aux[j].pt).squaredNorm();

        if (dist <= point.size * point.size) {
          if (point.response < kpts_aux[j].response) {
            is_repeated = true;
            break;
          }
        }
      }
    }

    if (is_repeated == false) {
      kpts.push_back(point);
    }
  }

  timing_.extrema = timer.elapsedMs();
}

/* ************************************************************************* */
void AKAZE::Do_Subpixel_Refinement(std::vector<AKAZEKeypoint>& kpts) {
  float Dx = 0.0, Dy = 0.0, ratio = 0.0;
  float Dxx = 0.0, Dyy = 0.0, Dxy = 0.0;
  int x = 0, y = 0;
  Eigen::Matrix2f A;
  A.setZero();
  Eigen::Vector2f b(0.0, 0.0);
  Eigen::Vector2f dst(0.0, 0.0);

  timer::Timer timer;
  for (size_t i = 0; i < kpts.size(); i++) {
    ratio = pow(2.f, kpts[i].octave);
    x = fRound(kpts[i].pt.x() / ratio);
    y = fRound(kpts[i].pt.y() / ratio);

    // Compute the gradient
    Dx = 0.5 * (evolution_[kpts[i].class_id].Ldet(y, x + 1) -
                evolution_[kpts[i].class_id].Ldet(y, x - 1));
    Dy = 0.5 * (evolution_[kpts[i].class_id].Ldet(y + 1, x) -
                evolution_[kpts[i].class_id].Ldet(y - 1, x));

    // Compute the Hessian
    Dxx = evolution_[kpts[i].class_id].Ldet(y, x + 1) +
          evolution_[kpts[i].class_id].Ldet(y, x - 1) -
          2.0 * evolution_[kpts[i].class_id].Ldet(y, x);

    Dyy = evolution_[kpts[i].class_id].Ldet(y + 1, x) +
          evolution_[kpts[i].class_id].Ldet(y - 1, x) -
          2.0 * evolution_[kpts[i].class_id].Ldet(y, x);

    Dxy = 0.25 * (evolution_[kpts[i].class_id].Ldet(y + 1, x + 1) +
                  evolution_[kpts[i].class_id].Ldet(y - 1, x - 1)) -
          0.25 * (evolution_[kpts[i].class_id].Ldet(y - 1, x + 1) +
                  evolution_[kpts[i].class_id].Ldet(y + 1, +x - 1));

    // Solve the linear system
    A(0, 0) = Dxx;
    A(1, 1) = Dyy;
    A(0, 1) = A(1, 0) = Dxy;
    b(0) = -Dx;
    b(1) = -Dy;

    dst = A.fullPivLu().solve(b);

    if (fabs(dst(0)) <= 1.0 && fabs(dst(1)) <= 1.0) {
      kpts[i].pt.x() = x + dst(0);
      kpts[i].pt.y() = y + dst(1);
      int power = powf(2, evolution_[kpts[i].class_id].octave);
      kpts[i].pt.x() *= power;
      kpts[i].pt.y() *= power;
      kpts[i].angle = 0.0;

      // In OpenCV the size of a keypoint its the diameter
      kpts[i].size *= 2.0;
    }
    // Delete the point since its not stable
    else {
      kpts.erase(kpts.begin() + i);
      i--;
    }
  }
  timing_.subpixel = timer.elapsedMs();
}

/* ************************************************************************* */
/**
 * @brief This method computes the set of descriptors through the nonlinear
 * scale space
 * @param kpts Vector of detected keypoints
 * @param desc Matrix to store the descriptors
*/
void AKAZE::Compute_Descriptors(std::vector<AKAZEKeypoint>& kpts,
                                AKAZEDescriptors& desc) {
  // Allocate memory for the matrix with the descriptors
  if (options_.descriptor < MLDB_UPRIGHT) {
    desc.float_descriptor.resize(kpts.size());
  } else {
    desc.binary_descriptor.resize(kpts.size());
    // We use the full length binary descriptor -> 486 bits
    int t;
    if (options_.descriptor_size == 0) {
      t = (6 + 36 + 120) * options_.descriptor_channels;
    } else {
      t = options_.descriptor_size;
    }

    for (int i = 0; i < desc.binary_descriptor.size(); i++) {
      desc.binary_descriptor[i].resize(t);
    }
  }

  timer::Timer timer;

  switch (options_.descriptor) {
    case SURF_UPRIGHT:  // Upright descriptors, not invariant to rotation
    {

#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)(kpts.size()); i++) {
        Get_SURF_Descriptor_Upright_64(kpts[i], desc.float_descriptor[i]);
      }
    } break;
    case SURF: {
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)(kpts.size()); i++) {
        Compute_Main_Orientation(kpts[i]);
        Get_SURF_Descriptor_64(kpts[i], desc.float_descriptor[i]);
      }
    } break;
    case MSURF_UPRIGHT:  // Upright descriptors, not invariant to rotation
    {
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)(kpts.size()); i++) {
        Get_MSURF_Upright_Descriptor_64(kpts[i], desc.float_descriptor[i]);
      }
    } break;
    case MSURF: {
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)(kpts.size()); i++) {
        Compute_Main_Orientation(kpts[i]);
        Get_MSURF_Descriptor_64(kpts[i], desc.float_descriptor[i]);
      }
    } break;
    case MLDB_UPRIGHT: {
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)(kpts.size()); i++) {
        Get_Upright_MLDB_Full_Descriptor(
            kpts[i], (unsigned char*)(desc.binary_descriptor[i].data()));
      }
    } break;
    case MLDB: {
#ifdef AKAZE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)(kpts.size()); i++) {
        Compute_Main_Orientation(kpts[i]);
        Get_MLDB_Full_Descriptor(
            kpts[i], (unsigned char*)(desc.binary_descriptor[i].data()));
      }
    } break;
  }

  timing_.descriptor = timer.elapsedMs();
}

/* ************************************************************************* */
void AKAZE::Compute_Main_Orientation(AKAZEKeypoint& kpt) const {
  int ix = 0, iy = 0, idx = 0, s = 0, level = 0;
  float xf = 0.0, yf = 0.0, gweight = 0.0, ratio = 0.0;
  float resX[109], resY[109], Ang[109];
  const int id[] = {6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6};

  // Variables for computing the dominant direction
  float sumX = 0.0, sumY = 0.0, max = 0.0, ang1 = 0.0, ang2 = 0.0;

  // Get the information from the keypoint
  level = kpt.class_id;
  ratio = (float)(1 << evolution_[level].octave);
  s = fRound(0.5 * kpt.size / ratio);
  xf = kpt.pt.x() / ratio;
  yf = kpt.pt.y() / ratio;

  // Calculate derivatives responses for points within radius of 6*scale
  for (int i = -6; i <= 6; ++i) {
    for (int j = -6; j <= 6; ++j) {
      if (i * i + j * j < 36) {
        iy = fRound(yf + j * s);
        ix = fRound(xf + i * s);

        gweight = gauss25[id[i + 6]][id[j + 6]];
        resX[idx] = gweight * evolution_[level].Lx(iy, ix);
        resY[idx] = gweight * evolution_[level].Ly(iy, ix);
        Ang[idx] = atan2(resY[idx], resX[idx]) * (M_PI / 180.0);
        ++idx;
      }
    }
  }

  // Loop slides pi/3 window around feature point
  for (ang1 = 0; ang1 < 2.0 * M_PI; ang1 += 0.15f) {
    ang2 = (ang1 + M_PI / 3.0f > 2.0 * M_PI ? ang1 - 5.0f * M_PI / 3.0f
                                            : ang1 + M_PI / 3.0f);
    sumX = sumY = 0.f;

    for (size_t k = 0; k < 109; ++k) {
      // Get angle from the x-axis of the sample point
      const float& ang = Ang[k];

      // Determine whether the point is within the window
      if (ang1 < ang2 && ang1 < ang && ang < ang2) {
        sumX += resX[k];
        sumY += resY[k];
      } else if (ang2 < ang1 && ((ang > 0 && ang < ang2) ||
                                 (ang > ang1 && ang < 2.0 * M_PI))) {
        sumX += resX[k];
        sumY += resY[k];
      }
    }

    // if the vector produced from this window is longer than all
    // previous vectors then this forms the new dominant direction
    if (sumX * sumX + sumY * sumY > max) {
      // store largest orientation
      max = sumX * sumX + sumY * sumY;
      kpt.angle = atan2(sumY, sumX) * (M_PI / 180.0);
    }
  }
}

/* ************************************************************************* */
void AKAZE::Get_SURF_Descriptor_Upright_64(const AKAZEKeypoint& kpt,
                                           Eigen::VectorXf& desc) const {
  desc.setZero(64);
  float dx = 0.0, dy = 0.0, mdx = 0.0, mdy = 0.0;
  float rx = 0.0, ry = 0.0, len = 0.0, xf = 0.0, yf = 0.0;
  float sample_x = 0.0, sample_y = 0.0;
  float fx = 0.0, fy = 0.0, ratio = 0.0, res1 = 0.0, res2 = 0.0, res3 = 0.0,
        res4 = 0.0;
  int x1 = 0, y1 = 0, x2 = 0, y2 = 0, sample_step = 0, pattern_size = 0,
      dcount = 0;
  int scale = 0, dsize = 0, level = 0;

  // Set the descriptor size and the sample and pattern sizes
  dsize = 64;
  sample_step = 5;
  pattern_size = 10;

  // Get the information from the keypoint
  ratio = (float)(1 << kpt.octave);
  scale = fRound(0.5 * kpt.size / ratio);
  level = kpt.class_id;
  yf = kpt.pt.y() / ratio;
  xf = kpt.pt.x() / ratio;

  // Calculate descriptor for this interest point
  for (int i = -pattern_size; i < pattern_size; i += sample_step) {
    for (int j = -pattern_size; j < pattern_size; j += sample_step) {

      dx = dy = mdx = mdy = 0.0;

      for (int k = i; k < i + sample_step; k++) {
        for (int l = j; l < j + sample_step; l++) {

          // Get the coordinates of the sample point on the rotated axis
          sample_y = yf + l * scale;
          sample_x = xf + k * scale;

          y1 = (int)(sample_y - .5);
          x1 = (int)(sample_x - .5);

          y2 = (int)(sample_y + .5);
          x2 = (int)(sample_x + .5);

          fx = sample_x - x1;
          fy = sample_y - y1;

          res1 = evolution_[level].Lx(y1, x1);
          res2 = evolution_[level].Lx(y1, x2);
          res3 = evolution_[level].Lx(y2, x1);
          res4 = evolution_[level].Lx(y2, x2);
          rx = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          res1 = evolution_[level].Ly(y1, x1);
          res2 = evolution_[level].Ly(y1, x2);
          res3 = evolution_[level].Ly(y2, x1);
          res4 = evolution_[level].Ly(y2, x2);
          ry = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          // Sum the derivatives to the cumulative descriptor
          dx += rx;
          dy += ry;
          mdx += fabs(rx);
          mdy += fabs(ry);
        }
      }

      // Add the values to the descriptor vector
      desc(dcount++) = dx;
      desc(dcount++) = dy;
      desc(dcount++) = mdx;
      desc(dcount++) = mdy;

      // Store the current length^2 of the vector
      len += dx * dx + dy * dy + mdx * mdx + mdy * mdy;
    }
  }

  // convert to unit vector
  len = sqrt(len);

  for (int i = 0; i < dsize; i++) {
    desc[i] /= len;
  }
}

/* ************************************************************************* */
void AKAZE::Get_SURF_Descriptor_64(const AKAZEKeypoint& kpt,
                                   Eigen::VectorXf& desc) const {
  desc.setZero(64);
  float dx = 0.0, dy = 0.0, mdx = 0.0, mdy = 0.0;
  float rx = 0.0, ry = 0.0, rrx = 0.0, rry = 0.0, len = 0.0, xf = 0.0, yf = 0.0;
  float sample_x = 0.0, sample_y = 0.0, co = 0.0, si = 0.0, angle = 0.0;
  float fx = 0.0, fy = 0.0, ratio = 0.0, res1 = 0.0, res2 = 0.0, res3 = 0.0,
        res4 = 0.0;
  int x1 = 0, y1 = 0, x2 = 0, y2 = 0, sample_step = 0, pattern_size = 0,
      dcount = 0;
  int scale = 0, dsize = 0, level = 0;

  // Set the descriptor size and the sample and pattern sizes
  dsize = 64;
  sample_step = 5;
  pattern_size = 10;

  // Get the information from the keypoint
  ratio = (float)(1 << kpt.octave);
  scale = fRound(0.5 * kpt.size / ratio);
  angle = kpt.angle;
  level = kpt.class_id;
  yf = kpt.pt.y() / ratio;
  xf = kpt.pt.x() / ratio;
  co = cos(angle);
  si = sin(angle);

  // Calculate descriptor for this interest point
  for (int i = -pattern_size; i < pattern_size; i += sample_step) {
    for (int j = -pattern_size; j < pattern_size; j += sample_step) {

      dx = dy = mdx = mdy = 0.0;

      for (int k = i; k < i + sample_step; k++) {
        for (int l = j; l < j + sample_step; l++) {

          // Get the coordinates of the sample point on the rotated axis
          sample_y = yf + (l * scale * co + k * scale * si);
          sample_x = xf + (-l * scale * si + k * scale * co);

          y1 = (int)(sample_y - .5);
          x1 = (int)(sample_x - .5);

          y2 = (int)(sample_y + .5);
          x2 = (int)(sample_x + .5);

          fx = sample_x - x1;
          fy = sample_y - y1;

          res1 = evolution_[level].Lx(y1, x1);
          res2 = evolution_[level].Lx(y1, x2);
          res3 = evolution_[level].Lx(y2, x1);
          res4 = evolution_[level].Lx(y2, x2);
          rx = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          res1 = evolution_[level].Ly(y1, x1);
          res2 = evolution_[level].Ly(y1, x2);
          res3 = evolution_[level].Ly(y2, x1);
          res4 = evolution_[level].Ly(y2, x2);
          ry = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          // Get the x and y derivatives on the rotated axis
          rry = rx * co + ry * si;
          rrx = -rx * si + ry * co;

          // Sum the derivatives to the cumulative descriptor
          dx += rrx;
          dy += rry;
          mdx += fabs(rrx);
          mdy += fabs(rry);
        }
      }

      // Add the values to the descriptor vector
      desc(dcount++) = dx;
      desc(dcount++) = dy;
      desc(dcount++) = mdx;
      desc(dcount++) = mdy;

      // Store the current length^2 of the vector
      len += dx * dx + dy * dy + mdx * mdx + mdy * mdy;
    }
  }

  // convert to unit vector
  len = sqrt(len);

  for (int i = 0; i < dsize; i++) desc[i] /= len;
}

/* ************************************************************************* */
void AKAZE::Get_MSURF_Upright_Descriptor_64(const AKAZEKeypoint& kpt,
                                            Eigen::VectorXf& desc) const {
  desc.setZero(64);
  float dx = 0.0, dy = 0.0, mdx = 0.0, mdy = 0.0, gauss_s1 = 0.0,
        gauss_s2 = 0.0;
  float rx = 0.0, ry = 0.0, len = 0.0, xf = 0.0, yf = 0.0, ys = 0.0, xs = 0.0;
  float sample_x = 0.0, sample_y = 0.0;
  int x1 = 0, y1 = 0, sample_step = 0, pattern_size = 0;
  int x2 = 0, y2 = 0, kx = 0, ky = 0, i = 0, j = 0, dcount = 0;
  float fx = 0.0, fy = 0.0, ratio = 0.0, res1 = 0.0, res2 = 0.0, res3 = 0.0,
        res4 = 0.0;
  int scale = 0, dsize = 0, level = 0;

  // Subregion centers for the 4x4 gaussian weighting
  float cx = -0.5, cy = 0.5;

  // Set the descriptor size and the sample and pattern sizes
  dsize = 64;
  sample_step = 5;
  pattern_size = 12;

  // Get the information from the keypoint
  ratio = (float)(1 << kpt.octave);
  scale = fRound(0.5 * kpt.size / ratio);
  level = kpt.class_id;
  yf = kpt.pt.y() / ratio;
  xf = kpt.pt.x() / ratio;

  i = -8;

  // Calculate descriptor for this interest point
  // Area of size 24 s x 24 s
  while (i < pattern_size) {
    j = -8;
    i = i - 4;

    cx += 1.0;
    cy = -0.5;

    while (j < pattern_size) {
      dx = dy = mdx = mdy = 0.0;
      cy += 1.0;
      j = j - 4;

      ky = i + sample_step;
      kx = j + sample_step;

      ys = yf + (ky * scale);
      xs = xf + (kx * scale);

      for (int k = i; k < i + 9; k++) {
        for (int l = j; l < j + 9; l++) {
          sample_y = k * scale + yf;
          sample_x = l * scale + xf;

          // Get the gaussian weighted x and y responses
          gauss_s1 = gaussian(xs - sample_x, ys - sample_y, 2.50 * scale);

          y1 = (int)(sample_y - .5);
          x1 = (int)(sample_x - .5);

          y2 = (int)(sample_y + .5);
          x2 = (int)(sample_x + .5);

          fx = sample_x - x1;
          fy = sample_y - y1;

          res1 = evolution_[level].Lx(y1, x1);
          res2 = evolution_[level].Lx(y1, x2);
          res3 = evolution_[level].Lx(y2, x1);
          res4 = evolution_[level].Lx(y2, x2);
          rx = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          res1 = evolution_[level].Ly(y1, x1);
          res2 = evolution_[level].Ly(y1, x2);
          res3 = evolution_[level].Ly(y2, x1);
          res4 = evolution_[level].Ly(y2, x2);
          ry = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          rx = gauss_s1 * rx;
          ry = gauss_s1 * ry;

          // Sum the derivatives to the cumulative descriptor
          dx += rx;
          dy += ry;
          mdx += fabs(rx);
          mdy += fabs(ry);
        }
      }

      // Add the values to the descriptor vector
      gauss_s2 = gaussian(cx - 2.0f, cy - 2.0f, 1.5f);

      desc(dcount++) = dx * gauss_s2;
      desc(dcount++) = dy * gauss_s2;
      desc(dcount++) = mdx * gauss_s2;
      desc(dcount++) = mdy * gauss_s2;

      len += (dx * dx + dy * dy + mdx * mdx + mdy * mdy) * gauss_s2 * gauss_s2;

      j += 9;
    }

    i += 9;
  }

  // convert to unit vector
  len = sqrt(len);

  for (int i = 0; i < dsize; i++) desc[i] /= len;
}

/* ************************************************************************* */
void AKAZE::Get_MSURF_Descriptor_64(const AKAZEKeypoint& kpt,
                                    Eigen::VectorXf& desc) const {
  desc.setZero(64);
  float dx = 0.0, dy = 0.0, mdx = 0.0, mdy = 0.0, gauss_s1 = 0.0,
        gauss_s2 = 0.0;
  float rx = 0.0, ry = 0.0, rrx = 0.0, rry = 0.0, len = 0.0, xf = 0.0, yf = 0.0,
        ys = 0.0, xs = 0.0;
  float sample_x = 0.0, sample_y = 0.0, co = 0.0, si = 0.0, angle = 0.0;
  float fx = 0.0, fy = 0.0, ratio = 0.0, res1 = 0.0, res2 = 0.0, res3 = 0.0,
        res4 = 0.0;
  int x1 = 0, y1 = 0, x2 = 0, y2 = 0, sample_step = 0, pattern_size = 0;
  int kx = 0, ky = 0, i = 0, j = 0, dcount = 0;
  int scale = 0, dsize = 0, level = 0;

  // Subregion centers for the 4x4 gaussian weighting
  float cx = -0.5, cy = 0.5;

  // Set the descriptor size and the sample and pattern sizes
  dsize = 64;
  sample_step = 5;
  pattern_size = 12;

  // Get the information from the keypoint
  ratio = (float)(1 << kpt.octave);
  scale = fRound(0.5 * kpt.size / ratio);
  angle = kpt.angle;
  level = kpt.class_id;
  yf = kpt.pt.y() / ratio;
  xf = kpt.pt.x() / ratio;
  co = cos(angle);
  si = sin(angle);

  i = -8;

  // Calculate descriptor for this interest point
  // Area of size 24 s x 24 s
  while (i < pattern_size) {
    j = -8;
    i = i - 4;

    cx += 1.0;
    cy = -0.5;

    while (j < pattern_size) {
      dx = dy = mdx = mdy = 0.0;
      cy += 1.0;
      j = j - 4;

      ky = i + sample_step;
      kx = j + sample_step;

      xs = xf + (-kx * scale * si + ky * scale * co);
      ys = yf + (kx * scale * co + ky * scale * si);

      for (int k = i; k < i + 9; ++k) {
        for (int l = j; l < j + 9; ++l) {
          // Get coords of sample point on the rotated axis
          sample_y = yf + (l * scale * co + k * scale * si);
          sample_x = xf + (-l * scale * si + k * scale * co);

          // Get the gaussian weighted x and y responses
          gauss_s1 = gaussian(xs - sample_x, ys - sample_y, 2.5 * scale);

          y1 = fRound(sample_y - .5);
          x1 = fRound(sample_x - .5);

          y2 = fRound(sample_y + .5);
          x2 = fRound(sample_x + .5);

          fx = sample_x - x1;
          fy = sample_y - y1;

          res1 = evolution_[level].Lx(y1, x1);
          res2 = evolution_[level].Lx(y1, x2);
          res3 = evolution_[level].Lx(y2, x1);
          res4 = evolution_[level].Lx(y2, x2);
          rx = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          res1 = evolution_[level].Ly(y1, x1);
          res2 = evolution_[level].Ly(y1, x2);
          res3 = evolution_[level].Ly(y2, x1);
          res4 = evolution_[level].Ly(y2, x2);
          ry = (1.0 - fx) * (1.0 - fy) * res1 + fx * (1.0 - fy) * res2 +
               (1.0 - fx) * fy * res3 + fx * fy * res4;

          // Get the x and y derivatives on the rotated axis
          rry = gauss_s1 * (rx * co + ry * si);
          rrx = gauss_s1 * (-rx * si + ry * co);

          // Sum the derivatives to the cumulative descriptor
          dx += rrx;
          dy += rry;
          mdx += fabs(rrx);
          mdy += fabs(rry);
        }
      }

      // Add the values to the descriptor vector
      gauss_s2 = gaussian(cx - 2.0f, cy - 2.0f, 1.5f);
      desc(dcount++) = dx * gauss_s2;
      desc(dcount++) = dy * gauss_s2;
      desc(dcount++) = mdx * gauss_s2;
      desc(dcount++) = mdy * gauss_s2;

      len += (dx * dx + dy * dy + mdx * mdx + mdy * mdy) * gauss_s2 * gauss_s2;

      j += 9;
    }

    i += 9;
  }

  // convert to unit vector
  len = sqrt(len);

  desc /= len;
}

/* ************************************************************************* */
void AKAZE::Get_MLDB_Full_Descriptor(const AKAZEKeypoint& kpt,
                                     unsigned char* desc) const {

  const int max_channels = 3;
  assert(options_.descriptor_channels <= max_channels);
  float values[16 * max_channels];
  const double size_mult[3] = {1, 2.0 / 3.0, 1.0 / 2.0};

  float ratio = (float)(1 << kpt.octave);
  float scale = (float)fRound(0.5f * kpt.size / ratio);
  float xf = kpt.pt.x() / ratio;
  float yf = kpt.pt.y() / ratio;
  float co = cos(kpt.angle);
  float si = sin(kpt.angle);
  int pattern_size = options_.descriptor_pattern_size;

  int dpos = 0;
  for (int lvl = 0; lvl < 3; lvl++) {
    int val_count = (lvl + 2) * (lvl + 2);
    int sample_step = static_cast<int>(ceil(pattern_size * size_mult[lvl]));
    MLDB_Fill_Values(values, sample_step, kpt.class_id, xf, yf, co, si, scale);
    MLDB_Binary_Comparisons(values, desc, val_count, dpos);
  }
}

/* ************************************************************************* */
void AKAZE::Get_Upright_MLDB_Full_Descriptor(const AKAZEKeypoint& kpt,
                                             unsigned char* desc) const {

  const int max_channels = 3;
  assert(options_.descriptor_channels <= max_channels);
  float values[16 * max_channels];
  const double size_mult[3] = {1, 2.0 / 3.0, 1.0 / 2.0};

  float ratio = (float)(1 << kpt.octave);
  float scale = (float)fRound(0.5f * kpt.size / ratio);
  float xf = kpt.pt.x() / ratio;
  float yf = kpt.pt.y() / ratio;
  int pattern_size = options_.descriptor_pattern_size;

  int dpos = 0;
  for (int lvl = 0; lvl < 3; lvl++) {
    int val_count = (lvl + 2) * (lvl + 2);
    int sample_step = static_cast<int>(ceil(pattern_size * size_mult[lvl]));
    MLDB_Fill_Upright_Values(values, sample_step, kpt.class_id, xf, yf, scale);
    MLDB_Binary_Comparisons(values, desc, val_count, dpos);
  }
}

/* ************************************************************************* */
void AKAZE::MLDB_Fill_Values(float* values, int sample_step, int level,
                             float xf, float yf, float co, float si,
                             float scale) const {

  int pattern_size = options_.descriptor_pattern_size;
  int nr_channels = options_.descriptor_channels;
  int valpos = 0;

  for (int i = -pattern_size; i < pattern_size; i += sample_step) {
    for (int j = -pattern_size; j < pattern_size; j += sample_step) {

      float di = 0.0, dx = 0.0, dy = 0.0;
      int nsamples = 0;

      for (int k = i; k < i + sample_step; k++) {
        for (int l = j; l < j + sample_step; l++) {

          float sample_y = yf + (l * co * scale + k * si * scale);
          float sample_x = xf + (-l * si * scale + k * co * scale);

          int y1 = fRound(sample_y);
          int x1 = fRound(sample_x);

          float ri = evolution_[level].Lt(y1, x1);
          di += ri;

          if (nr_channels > 1) {
            float rx = evolution_[level].Lx(y1, x1);
            float ry = evolution_[level].Ly(y1, x1);
            if (nr_channels == 2) {
              dx += sqrtf(rx * rx + ry * ry);
            } else {
              float rry = rx * co + ry * si;
              float rrx = -rx * si + ry * co;
              dx += rrx;
              dy += rry;
            }
          }
          nsamples++;
        }
      }

      di /= nsamples;
      dx /= nsamples;
      dy /= nsamples;

      values[valpos] = di;

      if (nr_channels > 1) values[valpos + 1] = dx;

      if (nr_channels > 2) values[valpos + 2] = dy;

      valpos += nr_channels;
    }
  }
}

/* ************************************************************************* */
void AKAZE::MLDB_Fill_Upright_Values(float* values, int sample_step, int level,
                                     float xf, float yf, float scale) const {

  int pattern_size = options_.descriptor_pattern_size;
  int nr_channels = options_.descriptor_channels;
  int valpos = 0;

  for (int i = -pattern_size; i < pattern_size; i += sample_step) {
    for (int j = -pattern_size; j < pattern_size; j += sample_step) {

      float di = 0.0, dx = 0.0, dy = 0.0;
      int nsamples = 0;

      for (int k = i; k < i + sample_step; k++) {
        for (int l = j; l < j + sample_step; l++) {

          float sample_y = yf + l * scale;
          float sample_x = xf + k * scale;

          int y1 = fRound(sample_y);
          int x1 = fRound(sample_x);

          float ri = evolution_[level].Lt(y1, x1);
          di += ri;

          if (nr_channels > 1) {
            float rx = evolution_[level].Lx(y1, x1);
            float ry = evolution_[level].Ly(y1, x1);
            if (nr_channels == 2) {
              dx += sqrtf(rx * rx + ry * ry);
            } else {
              dx += rx;
              dy += ry;
            }
          }
          nsamples++;
        }
      }

      di /= nsamples;
      dx /= nsamples;
      dy /= nsamples;

      values[valpos] = di;

      if (nr_channels > 1) values[valpos + 1] = dx;

      if (nr_channels > 2) values[valpos + 2] = dy;

      valpos += nr_channels;
    }
  }
}

/* ************************************************************************* */
void AKAZE::MLDB_Binary_Comparisons(float* values, unsigned char* desc,
                                    int count, int& dpos) const {
#define TOGGLE_FLT(x) ((x) ^ ((int)(x) < 0 ? 0x7fffffff : 0))

  int nr_channels = options_.descriptor_channels;
  int* ivalues = (int*)values;
  for (int i = 0; i < count * nr_channels; i++) {
    ivalues[i] = TOGGLE_FLT(ivalues[i]);
  }
#undef TOGGLE_FLT

  for (int pos = 0; pos < nr_channels; pos++) {
    for (int i = 0; i < count; i++) {
      int ival = ivalues[nr_channels * i + pos];
      for (int j = i + 1; j < count; j++) {
        int res = ival > ivalues[nr_channels * j + pos];
        desc[dpos >> 3] |= (res << (dpos & 7));
        dpos++;
      }
    }
  }
}

/* ************************************************************************* */
void AKAZE::Get_MLDB_Descriptor_Subset(const AKAZEKeypoint& kpt,
                                       unsigned char* desc) {
  float di = 0.f, dx = 0.f, dy = 0.f;
  float rx = 0.f, ry = 0.f;
  float sample_x = 0.f, sample_y = 0.f;
  int x1 = 0, y1 = 0;

  // Get the information from the keypoint
  float ratio = (float)(1 << kpt.octave);
  int scale = fRound(0.5 * kpt.size / ratio);
  float angle = kpt.angle;
  float level = kpt.class_id;
  float yf = kpt.pt.y() / ratio;
  float xf = kpt.pt.x() / ratio;
  float co = cos(angle);
  float si = sin(angle);

  // Allocate memory for the matrix of values
  Eigen::VectorXf values((4 + 9 + 16) * options_.descriptor_channels);

  // Sample everything, but only do the comparisons
  std::vector<int> steps(3);
  steps.at(0) = options_.descriptor_pattern_size;
  steps.at(1) = ceil(2.f * options_.descriptor_pattern_size / 3.f);
  steps.at(2) = options_.descriptor_pattern_size / 2;

  for (int i = 0; i < descriptorSamples_.rows(); i++) {
    const Eigen::RowVectorXi& coords = descriptorSamples_.row(i);
    int sample_step = steps.at(coords(0));
    di = 0.0f;
    dx = 0.0f;
    dy = 0.0f;

    for (int k = coords(1); k < coords(1) + sample_step; k++) {
      for (int l = coords(2); l < coords(2) + sample_step; l++) {

        // Get the coordinates of the sample point
        sample_y = yf + (l * scale * co + k * scale * si);
        sample_x = xf + (-l * scale * si + k * scale * co);

        y1 = fRound(sample_y);
        x1 = fRound(sample_x);

        di += evolution_[level].Lt(y1, x1);

        if (options_.descriptor_channels > 1) {
          rx = evolution_[level].Lx(y1, x1);
          ry = evolution_[level].Ly(y1, x1);

          if (options_.descriptor_channels == 2) {
           dx += sqrtf(rx * rx + ry * ry);
          } else if (options_.descriptor_channels == 3) {
            // Get the x and y derivatives on the rotated axis
            dx += rx * co + ry * si;
            dy += -rx * si + ry * co;
          }
        }
      }
    }

    values(options_.descriptor_channels * i) = di;

    if (options_.descriptor_channels == 2) {
      values(options_.descriptor_channels * i + 1) = dx;
    } else if (options_.descriptor_channels == 3) {
      values(options_.descriptor_channels * i + 1) = dx;
      values(options_.descriptor_channels * i + 2) = dy;
    }
  }

  // Do the comparisons
  const float* vals = values.data();
  const int* comps = descriptorBits_.data();

  for (int i = 0; i < descriptorBits_.rows(); i++) {
    if (vals[comps[2 * i]] > vals[comps[2 * i + 1]]) {
      desc[i / 8] |= (1 << (i % 8));
    }
  }
}

/* ************************************************************************* */
void AKAZE::Get_Upright_MLDB_Descriptor_Subset(const AKAZEKeypoint& kpt,
                                               unsigned char* desc) {

  float di = 0.0f, dx = 0.0f, dy = 0.0f;
  float rx = 0.0f, ry = 0.0f;
  float sample_x = 0.0f, sample_y = 0.0f;
  int x1 = 0, y1 = 0;

  // Get the information from the keypoint
  float ratio = (float)(1 << kpt.octave);
  int scale = fRound(0.5 * kpt.size / ratio);
  float level = kpt.class_id;
  float yf = kpt.pt.y() / ratio;
  float xf = kpt.pt.x() / ratio;

  // Allocate memory for the matrix of values
  Eigen::VectorXf values((4 + 9 + 16) * options_.descriptor_channels);

  std::vector<int> steps(3);
  steps.at(0) = options_.descriptor_pattern_size;
  steps.at(1) = ceil(2.f * options_.descriptor_pattern_size / 3.f);
  steps.at(2) = options_.descriptor_pattern_size / 2;

  for (int i = 0; i < descriptorSamples_.rows(); i++) {
    const Eigen::RowVectorXi coords = descriptorSamples_.row(i);
    int sample_step = steps.at(coords[0]);
    di = 0.0f, dx = 0.0f, dy = 0.0f;

    for (int k = coords[1]; k < coords[1] + sample_step; k++) {
      for (int l = coords[2]; l < coords[2] + sample_step; l++) {

        // Get the coordinates of the sample point
        sample_y = yf + l * scale;
        sample_x = xf + k * scale;

        y1 = fRound(sample_y);
        x1 = fRound(sample_x);
        di += evolution_[level].Lt(y1, x1);

        if (options_.descriptor_channels > 1) {
          rx = evolution_[level].Lx(y1, x1);
          ry = evolution_[level].Ly(y1, x1);

          if (options_.descriptor_channels == 2) {
            dx += sqrtf(rx * rx + ry * ry);
          } else if (options_.descriptor_channels == 3) {
            dx += rx;
            dy += ry;
          }
        }
      }
    }

    values(options_.descriptor_channels * i) = di;

    if (options_.descriptor_channels == 2) {
      values(options_.descriptor_channels * i + 1) = dx;
    } else if (options_.descriptor_channels == 3) {
      values(options_.descriptor_channels * i + 1) = dx;
      values(options_.descriptor_channels * i + 2) = dy;
    }
  }

  // Do the comparisons
  const float* vals = values.data();
  const int* comps = descriptorBits_.data();
  for (int i = 0; i < descriptorBits_.rows(); i++) {
    if (vals[comps[2 * i]] > vals[comps[2 * i + 1]]) {
      desc[i / 8] |= (1 << (i % 8));
    }
  }
}

/* ************************************************************************* */

// This method saves the scale space into jpg images
void AKAZE::Save_Scale_Space() {
  std::string outputFile;
  cimg_library::CImg<float> img_aux;
  for (size_t i = 0; i < evolution_.size(); i++) {
    ConvertEigenToCImg(evolution_[i].Lt, img_aux);
    outputFile = "../output/evolution_" + to_formatted_string(i, 2) + ".jpg";
    img_aux.save(outputFile.c_str());
  }
}

/* ************************************************************************* */

// This method saves the feature detector responses of the nonlinear scale
// space into jpg images
void AKAZE::Save_Detector_Responses() {
  cimg_library::CImg<float> img_aux;
  std::string outputFile;
  float ttime = 0.0;
  int nimgs = 0;

  for (size_t i = 0; i < evolution_.size(); i++) {
    ttime = evolution_[i + 1].etime - evolution_[i].etime;
    if (ttime > 0) {
      ConvertEigenToCImg(evolution_[i].Ldet, img_aux);
      outputFile =
          "../output/images/detector_" + to_formatted_string(nimgs, 2) + ".jpg";
      img_aux.save(outputFile.c_str());
      nimgs++;
    }
  }
}

/* ************************************************************************* */
void AKAZE::Show_Computation_Times() const {
  std::cout << "(*) Time Scale Space: " << timing_.scale << std::endl;
  std::cout << "(*) Time Detector: " << timing_.detector << std::endl;
  std::cout << "   - Time Derivatives: " << timing_.derivatives << std::endl;
  std::cout << "   - Time Extrema: " << timing_.extrema << std::endl;
  std::cout << "   - Time Subpixel: " << timing_.subpixel << std::endl;
  std::cout << "(*) Time Descriptor: " << timing_.descriptor << std::endl;
  std::cout << std::endl;
}

/* ************************************************************************* */
void libAKAZE::generateDescriptorSubsample(Eigen::MatrixXi& sampleList,
                                           Eigen::MatrixXi& comparisons, int nbits,
                                           int pattern_size, int nchannels) {
  int ssz = 0;
  for (int i = 0; i < 3; i++) {
    int gz = (i + 2) * (i + 2);
    ssz += gz * (gz - 1) / 2;
  }
  ssz *= nchannels;

  assert(nbits <= ssz &&
         "descriptor size can't be bigger than full descriptor");

  // Since the full descriptor is usually under 10k elements, we pick
  // the selection from the full matrix.  We take as many samples per
  // pick as the number of channels. For every pick, we
  // take the two samples involved and put them in the sampling list
  Eigen::MatrixXi fullM(ssz / nchannels, 5);
  for (size_t i = 0, c = 0; i < 3; i++) {
    const int gdiv = i + 2;  // grid divisions, per row
    const int gsz = gdiv * gdiv;
    const int psz = ceil(2. * pattern_size / (float)gdiv);

    for (int j = 0; j < gsz; j++) {
      for (int k = j + 1; k < gsz; k++, c++) {
        fullM(c, 0) = i;
        fullM(c, 1) = psz * (j % gdiv) - pattern_size;
        fullM(c, 2) = psz * (j / gdiv) - pattern_size;
        fullM(c, 3) = psz * (k % gdiv) - pattern_size;
        fullM(c, 4) = psz * (k / gdiv) - pattern_size;
      }
    }
  }

  srand(1024);
  Eigen::MatrixXi comps(
      nchannels * static_cast<int>(ceil(nbits / (float)nchannels)), 2);
  comps.setConstant(1000);

  // Select some samples. A sample includes all channels
  int count = 0;
  size_t npicks = ceil(nbits / (float)nchannels);
  Eigen::MatrixXi samples(29, 3);
  Eigen::MatrixXi fullcopy = fullM;
  samples.setConstant(-1);

  for (size_t i = 0; i < npicks; i++) {
    size_t k = rand() % (fullM.rows() - i);
    if (i < 6) {
      // Force use of the coarser grid values and comparisons
      k = i;
    }

    bool n = true;

    for (int j = 0; j < count; j++) {
      if (samples(j, 0) == fullcopy(k, 0) && samples(j, 1) == fullcopy(k, 1) &&
          samples(j, 2) == fullcopy(k, 2)) {
        n = false;
        comps(i * nchannels, 0) = nchannels * j;
        comps(i * nchannels + 1, 0) = nchannels * j + 1;
        comps(i * nchannels + 2, 0) = nchannels * j + 2;
        break;
      }
    }

    if (n) {
      samples(count, 0) = fullcopy(k, 0);
      samples(count, 1) = fullcopy(k, 1);
      samples(count, 2) = fullcopy(k, 2);
      comps(i * nchannels, 0) = nchannels * count;
      comps(i * nchannels + 1, 0) = nchannels * count + 1;
      comps(i * nchannels + 2, 0) = nchannels * count + 2;
      count++;
    }

    n = true;
    for (int j = 0; j < count; j++) {
      if (samples(j, 0) == fullcopy(k, 0) && samples(j, 1) == fullcopy(k, 3) &&
          samples(j, 2) == fullcopy(k, 4)) {
        n = false;
        comps(i * nchannels, 1) = nchannels * j;
        comps(i * nchannels + 1, 1) = nchannels * j + 1;
        comps(i * nchannels + 2, 1) = nchannels * j + 2;
        break;
      }
    }

    if (n) {
      samples(count, 0) = fullcopy(k, 0);
      samples(count, 1) = fullcopy(k, 3);
      samples(count, 2) = fullcopy(k, 4);
      comps(i * nchannels, 1) = nchannels * count;
      comps(i * nchannels + 1, 1) = nchannels * count + 1;
      comps(i * nchannels + 2, 1) = nchannels * count + 2;
      count++;
    }

    fullcopy.row(k) = fullcopy.row(fullcopy.rows() - i - 1);
  }

  sampleList = samples.topRows(count);
  comparisons = comps.topRows(nbits);
}

/* ************************************************************************* */
void libAKAZE::check_descriptor_limits(int& x, int& y, int width, int height) {

  if (x < 0) x = 0;

  if (y < 0) y = 0;

  if (x > width - 1) x = width - 1;

  if (y > height - 1) y = height - 1;
}
