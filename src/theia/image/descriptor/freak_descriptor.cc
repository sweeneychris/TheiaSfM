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
//
//  License of original FREAK code:
//
// Copyright (C) 2011-2012  Signal processing laboratory 2, EPFL,
// Kirell Benzi (kirell.benzi@epfl.ch),
// Raphael Ortiz (raphael.ortiz@a3.epfl.ch)
// Alexandre Alahi (alexandre.alahi@epfl.ch)
// and Pierre Vandergheynst (pierre.vandergheynst@epfl.ch)
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote
//     products derived from this software without specific prior written
//     permission.
//
//  This software is provided by the copyright holders and contributors "as is"
//  and any express or implied warranties, including, but not limited to, the
//  implied warranties of merchantability and fitness for a particular purpose
//  are disclaimed.  In no event shall the Intel Corporation or contributors be
//  liable for any direct, indirect, incidental, special, exemplary, or
//  consequential damages (including, but not limited to, procurement of
//  substitute goods or services; loss of use, data, or profits; or business
//  interruption) however caused and on any theory of liability, whether in
//  contract, strict liability, or tort (including negligence or otherwise)
//  arising in any way out of the use of this software, even if advised of the
//  possibility of such damage.

#include "theia/image/descriptor/freak_descriptor.h"

#define _USE_MATH_DEFINES
#include <glog/logging.h>
#include <stdint.h>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <vector>

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/image.h"
#include "theia/image/keypoint_detector/brisk_detector.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/util/util.h"

// It should be noted that this implementation is heavily derived from the
// implementation provided by the authors online at
// http://www.ivpe.com/freak.htm It has been modified significantly to reflect
// the use case for this library (e.g. we assume the default patterns given by
// the authors, rather than train them ourselves).
namespace theia {
namespace {

static const double kLog2 = 0.693147180559945;
static const int kNumOrientation = 256;
static const int kNumPoints = 43;
static const int kSmallestKeypointSize = 7;
// duplicate of Freak member var...
static const int kNumPairs = 512;
static const int freak_def_pairs[kNumPairs] = {
  404, 431, 818, 511, 181, 52, 311, 874, 774, 543, 719, 230, 417, 205, 11, 560,
  149, 265, 39, 306, 165, 857, 250, 8, 61, 15, 55, 717, 44, 412, 592, 134, 761,
  695, 660, 782, 625, 487, 549, 516, 271, 665, 762, 392, 178, 796, 773, 31, 672,
  845, 548, 794, 677, 654, 241, 831, 225, 238, 849, 83, 691, 484, 826, 707, 122,
  517, 583, 731, 328, 339, 571, 475, 394, 472, 580, 381, 137, 93, 380, 327, 619,
  729, 808, 218, 213, 459, 141, 806, 341, 95, 382, 568, 124, 750, 193, 749, 706,
  843, 79, 199, 317, 329, 768, 198, 100, 466, 613, 78, 562, 783, 689, 136, 838,
  94, 142, 164, 679, 219, 419, 366, 418, 423, 77, 89, 523, 259, 683, 312, 555,
  20, 470, 684, 123, 458, 453, 833, 72, 113, 253, 108, 313, 25, 153, 648, 411,
  607, 618, 128, 305, 232, 301, 84, 56, 264, 371, 46, 407, 360, 38, 99, 176,
  710, 114, 578, 66, 372, 653, 129, 359, 424, 159, 821, 10, 323, 393, 5, 340,
  891, 9, 790, 47, 0, 175, 346, 236, 26, 172, 147, 574, 561, 32, 294, 429, 724,
  755, 398, 787, 288, 299, 769, 565, 767, 722, 757, 224, 465, 723, 498, 467,
  235, 127, 802, 446, 233, 544, 482, 800, 318, 16, 532, 801, 441, 554, 173, 60,
  530, 713, 469, 30, 212, 630, 899, 170, 266, 799, 88, 49, 512, 399, 23, 500,
  107, 524, 90, 194, 143, 135, 192, 206, 345, 148, 71, 119, 101, 563, 870, 158,
  254, 214, 276, 464, 332, 725, 188, 385, 24, 476, 40, 231, 620, 171, 258, 67,
  109, 844, 244, 187, 388, 701, 690, 50, 7, 850, 479, 48, 522, 22, 154, 12, 659,
  736, 655, 577, 737, 830, 811, 174, 21, 237, 335, 353, 234, 53, 270, 62, 182,
  45, 177, 245, 812, 673, 355, 556, 612, 166, 204, 54, 248, 365, 226, 242, 452,
  700, 685, 573, 14, 842, 481, 468, 781, 564, 416, 179, 405, 35, 819, 608, 624,
  367, 98, 643, 448, 2, 460, 676, 440, 240, 130, 146, 184, 185, 430, 65, 807,
  377, 82, 121, 708, 239, 310, 138, 596, 730, 575, 477, 851, 797, 247, 27, 85,
  586, 307, 779, 326, 494, 856, 324, 827, 96, 748, 13, 397, 125, 688, 702, 92,
  293, 716, 277, 140, 112, 4, 80, 855, 839, 1, 413, 347, 584, 493, 289, 696, 19,
  751, 379, 76, 73, 115, 6, 590, 183, 734, 197, 483, 217, 344, 330, 400, 186,
  243, 587, 220, 780, 200, 793, 246, 824, 41, 735, 579, 81, 703, 322, 760, 720,
  139, 480, 490, 91, 814, 813, 163, 152, 488, 763, 263, 425, 410, 576, 120, 319,
  668, 150, 160, 302, 491, 515, 260, 145, 428, 97, 251, 395, 272, 252, 18, 106,
  358, 854, 485, 144, 550, 131, 133, 378, 68, 102, 104, 58, 361, 275, 209, 697,
  582, 338, 742, 589, 325, 408, 229, 28, 304, 191, 189, 110, 126, 486, 211, 547,
  533, 70, 215, 670, 249, 36, 581, 389, 605, 331, 518, 442, 822
};

}  // namespace

// Initializes the sampling patterns and local variables.
bool FreakDescriptorExtractor::Initialize() {
  // If it has already been initialized return.
  if (!pattern_lookup_.empty()) return true;

  pattern_lookup_.resize(kNumScales_ * kNumOrientation * kNumPoints);
  // 2 ^  ((num_octaves_-1) /nbScales)
  double scaleStep =
      std::pow(2.0, static_cast<double>(num_octaves_) / kNumScales_);
  double scalingFactor, alpha, beta, theta = 0;

  // pattern definition, radius normalized to 1.0 (outer point
  // position+sigma=1.0) number of points on each concentric circle (from outer
  // to inner)
  const int n[8] = { 6, 6, 6, 6, 6, 6, 6, 1 };
  // bigger radius
  const double bigR(2.0 / 3.0);
  // smaller radius
  const double smallR(2.0 / 24.0);
  // define spaces between concentric circles (from center to outer:
  // 1,2,3,4,5,6)
  const double unitSpace((bigR - smallR) / 21.0);
  // radii of the concentric cirles (from outer to inner)
  const double radius[8] = { bigR, bigR - 6 * unitSpace, bigR - 11 * unitSpace,
                             bigR - 15 * unitSpace, bigR - 18 * unitSpace,
                             bigR - 20 * unitSpace, smallR, 0.0 };
  // sigma of pattern points (each group of 6 points on a concentric cirle has
  // the same sigma)
  const double sigma[8] = { radius[0] / 2.0, radius[1] / 2.0, radius[2] / 2.0,
                            radius[3] / 2.0, radius[4] / 2.0, radius[5] / 2.0,
                            radius[6] / 2.0, radius[6] / 2.0 };
  // fill the lookup table
  for (int scaleIdx = 0; scaleIdx < kNumScales_; ++scaleIdx) {
    // proper initialization
    pattern_sizes_[scaleIdx] = 0;
    // scale of the pattern, scaleStep ^ scaleIdx
    scalingFactor = std::pow(scaleStep, scaleIdx);

    for (int orientationIdx = 0; orientationIdx < kNumOrientation;
         ++orientationIdx) {
      // orientation of the pattern
      theta = static_cast<double>(orientationIdx) * 2 * M_PI /
              static_cast<double>(kNumOrientation);
      int pointIdx = 0;

      for (size_t i = 0; i < 8; ++i) {
        for (int k = 0; k < n[i]; ++k) {
          // orientation offset so that groups of points on each circles are
          // staggered
          beta = M_PI / n[i] * (i % 2);
          alpha = k * 2.0 * M_PI / static_cast<double>(n[i]) + beta + theta;

          // add the point to the look-up table
          PatternPoint& point =
              pattern_lookup_[scaleIdx * kNumOrientation * kNumPoints +
                              orientationIdx * kNumPoints + pointIdx];
          point.x = static_cast<float>(radius[i] * cos(alpha) * scalingFactor *
                                       pattern_scale_);
          point.y = static_cast<float>(radius[i] * sin(alpha) * scalingFactor *
                                       pattern_scale_);
          point.sigma =
              static_cast<float>(sigma[i] * scalingFactor * pattern_scale_);

          // adapt the sizeList if necessary
          const int sizeMax = static_cast<int>(ceil(
              (radius[i] + sigma[i]) * scalingFactor * pattern_scale_)) + 1;
          if (pattern_sizes_[scaleIdx] < sizeMax)
            pattern_sizes_[scaleIdx] = sizeMax;

          ++pointIdx;
        }
      }
    }
  }

  // build the list of orientation pairs
  orientation_pairs_[0].i = 0;
  orientation_pairs_[0].j = 3;
  orientation_pairs_[1].i = 1;
  orientation_pairs_[1].j = 4;
  orientation_pairs_[2].i = 2;
  orientation_pairs_[2].j = 5;
  orientation_pairs_[3].i = 0;
  orientation_pairs_[3].j = 2;
  orientation_pairs_[4].i = 1;
  orientation_pairs_[4].j = 3;
  orientation_pairs_[5].i = 2;
  orientation_pairs_[5].j = 4;
  orientation_pairs_[6].i = 3;
  orientation_pairs_[6].j = 5;
  orientation_pairs_[7].i = 4;
  orientation_pairs_[7].j = 0;
  orientation_pairs_[8].i = 5;
  orientation_pairs_[8].j = 1;

  orientation_pairs_[9].i = 6;
  orientation_pairs_[9].j = 9;
  orientation_pairs_[10].i = 7;
  orientation_pairs_[10].j = 10;
  orientation_pairs_[11].i = 8;
  orientation_pairs_[11].j = 11;
  orientation_pairs_[12].i = 6;
  orientation_pairs_[12].j = 8;
  orientation_pairs_[13].i = 7;
  orientation_pairs_[13].j = 9;
  orientation_pairs_[14].i = 8;
  orientation_pairs_[14].j = 10;
  orientation_pairs_[15].i = 9;
  orientation_pairs_[15].j = 11;
  orientation_pairs_[16].i = 10;
  orientation_pairs_[16].j = 6;
  orientation_pairs_[17].i = 11;
  orientation_pairs_[17].j = 7;

  orientation_pairs_[18].i = 12;
  orientation_pairs_[18].j = 15;
  orientation_pairs_[19].i = 13;
  orientation_pairs_[19].j = 16;
  orientation_pairs_[20].i = 14;
  orientation_pairs_[20].j = 17;
  orientation_pairs_[21].i = 12;
  orientation_pairs_[21].j = 14;
  orientation_pairs_[22].i = 13;
  orientation_pairs_[22].j = 15;
  orientation_pairs_[23].i = 14;
  orientation_pairs_[23].j = 16;
  orientation_pairs_[24].i = 15;
  orientation_pairs_[24].j = 17;
  orientation_pairs_[25].i = 16;
  orientation_pairs_[25].j = 12;
  orientation_pairs_[26].i = 17;
  orientation_pairs_[26].j = 13;

  orientation_pairs_[27].i = 18;
  orientation_pairs_[27].j = 21;
  orientation_pairs_[28].i = 19;
  orientation_pairs_[28].j = 22;
  orientation_pairs_[29].i = 20;
  orientation_pairs_[29].j = 23;
  orientation_pairs_[30].i = 18;
  orientation_pairs_[30].j = 20;
  orientation_pairs_[31].i = 19;
  orientation_pairs_[31].j = 21;
  orientation_pairs_[32].i = 20;
  orientation_pairs_[32].j = 22;
  orientation_pairs_[33].i = 21;
  orientation_pairs_[33].j = 23;
  orientation_pairs_[34].i = 22;
  orientation_pairs_[34].j = 18;
  orientation_pairs_[35].i = 23;
  orientation_pairs_[35].j = 19;

  orientation_pairs_[36].i = 24;
  orientation_pairs_[36].j = 27;
  orientation_pairs_[37].i = 25;
  orientation_pairs_[37].j = 28;
  orientation_pairs_[38].i = 26;
  orientation_pairs_[38].j = 29;
  orientation_pairs_[39].i = 30;
  orientation_pairs_[39].j = 33;
  orientation_pairs_[40].i = 31;
  orientation_pairs_[40].j = 34;
  orientation_pairs_[41].i = 32;
  orientation_pairs_[41].j = 35;
  orientation_pairs_[42].i = 36;
  orientation_pairs_[42].j = 39;
  orientation_pairs_[43].i = 37;
  orientation_pairs_[43].j = 40;
  orientation_pairs_[44].i = 38;
  orientation_pairs_[44].j = 41;

  for (unsigned m = kNumOrientationPairs_; m--;) {
    const float dx = pattern_lookup_[orientation_pairs_[m].i].x -
                     pattern_lookup_[orientation_pairs_[m].j].x;
    const float dy = pattern_lookup_[orientation_pairs_[m].i].y -
                     pattern_lookup_[orientation_pairs_[m].j].y;
    const float norm_sq = (dx * dx + dy * dy);
    orientation_pairs_[m].weight_dx =
        static_cast<int>((dx / (norm_sq)) * 4096.0 + 0.5);
    orientation_pairs_[m].weight_dy =
        static_cast<int>((dy / (norm_sq)) * 4096.0 + 0.5);
  }

  // build the list of description pairs
  std::vector<DescriptionPair> allPairs;
  for (unsigned int i = 1; i < (unsigned int) kNumPoints; ++i) {
    // (generate all the pairs)
    for (unsigned int j = 0; (unsigned int) j < i; ++j) {
      DescriptionPair pair = {(uchar) i, (uchar) j };
      allPairs.push_back(pair);
    }
  }
  for (int i = 0; i < kNumPairs_; ++i)
    description_pairs_[i] = allPairs[freak_def_pairs[i]];
  return true;
}

// Computes a descriptor at a single keypoint.
bool FreakDescriptorExtractor::ComputeDescriptor(
    const FloatImage& image,
    const Keypoint& keypoint,
    BinaryVectorX* descriptor) {
  std::vector<Keypoint> keypoints;
  keypoints.push_back(keypoint);
  std::vector<BinaryVectorX> descriptors;
  bool success =
      ComputeDescriptors(image, &keypoints, &descriptors);
  if (!success || keypoints.size() == 0) {
    return false;
  }

  *descriptor = descriptors[0];
  return true;
}

// Compute multiple descriptors for keypoints from a single image.
bool FreakDescriptorExtractor::ComputeDescriptors(
    const FloatImage& image,
    std::vector<Keypoint>* keypoints,
    std::vector<BinaryVectorX>* descriptors) {
  Image<uchar> uchar_image(image.AsGrayscaleImage());
  Image<int> img_integral;
  uchar_image.Integrate(&img_integral);

  // used to save pattern scale index corresponding to each keypoints
  std::vector<int> kp_scale_idx(keypoints->size());
  const float size_cst =
      static_cast<float>(kNumScales_ / (kLog2 * num_octaves_));
  uchar points_value[kNumPoints];
  int theta_idx = 0;
  int direction0;
  int direction1;

  // compute the scale index corresponding to the keypoint size and remove
  // keypoints close to the border.
  descriptors->reserve(keypoints->size());
  if (scale_invariant_) {
    for (size_t k = keypoints->size(); k--;) {
      // Is k non-zero? If so, decrement it and continue.
      kp_scale_idx[k] = std::max(
          static_cast<int>(std::log((*keypoints)[k].scale() /
                                    kSmallestKeypointSize) * size_cst + 0.5),
          0);
      if (kp_scale_idx[k] >= kNumScales_) {
        kp_scale_idx[k] = kNumScales_ - 1;
      }

      // Check if the description at this specific position and scale fits
      // inside the image.
      if ((*keypoints)[k].x() <= pattern_sizes_[kp_scale_idx[k]] ||
          (*keypoints)[k].y() <= pattern_sizes_[kp_scale_idx[k]] ||
          (*keypoints)[k].x() >=
              image.Cols() - pattern_sizes_[kp_scale_idx[k]] ||
          (*keypoints)[k].y() >=
              image.Rows() - pattern_sizes_[kp_scale_idx[k]]) {
        keypoints->erase(keypoints->begin() + k);
        kp_scale_idx.erase(kp_scale_idx.begin() + k);
      }
    }
  } else {
    // Equivalent to the formula when the scale is normalized with a constant
    // size of (*keypoints)[k].size=3*SMALLEST_KP_SIZE.
    int scIdx = std::max(static_cast<int>(1.0986122886681 * size_cst + 0.5), 0);
    if (scIdx >= kNumScales_) {
      scIdx = kNumScales_ - 1;
    }
    for (size_t k = keypoints->size(); k--;) {
      kp_scale_idx[k] = scIdx;
      if ((*keypoints)[k].x() <= pattern_sizes_[kp_scale_idx[k]] ||
          (*keypoints)[k].y() <= pattern_sizes_[kp_scale_idx[k]] ||
          (*keypoints)[k].x() >=
              image.Cols() - pattern_sizes_[kp_scale_idx[k]] ||
          (*keypoints)[k].y() >=
              image.Rows() - pattern_sizes_[kp_scale_idx[k]]) {
        keypoints->erase(keypoints->begin() + k);
        kp_scale_idx.erase(kp_scale_idx.begin() + k);
      }
    }
  }

  // Estimate orientations, extract descriptors, extract the best comparisons
  // only.
  for (size_t k = keypoints->size(); k--;) {
    BinaryVectorX freak_descriptor(512 / (8 * sizeof(uint8_t)));
    // estimate orientation (gradient)
    if (!rotation_invariant_) {
      // assign 0° to all keypoints
      theta_idx = 0;
    } else {
      // get the points intensity value in the un-rotated pattern
      for (int i = kNumPoints; i--;) {
        points_value[i] = MeanIntensity(
            uchar_image, img_integral, keypoints->at(k).x(),
            keypoints->at(k).y(), kp_scale_idx[k], 0, i);
      }
      direction0 = 0;
      direction1 = 0;
      for (int m = 45; m--;) {
        // iterate through the orientation pairs
        const int delta = points_value[orientation_pairs_[m].i] -
                          points_value[orientation_pairs_[m].j];
        direction0 += delta * (orientation_pairs_[m].weight_dx) / 2048;
        direction1 += delta * (orientation_pairs_[m].weight_dy) / 2048;
      }
      // estimate orientation
      const float orientation = static_cast<float>(
          atan2(static_cast<float>(direction1),
                static_cast<float>(direction0) * (180.0 / M_PI)));
      theta_idx =
          static_cast<int>(kNumOrientation * orientation * (1 / 360.0) + 0.5);
      if (theta_idx < 0) theta_idx += kNumOrientation;

      if (theta_idx >= kNumOrientation) theta_idx -= kNumOrientation;
    }
    // extract descriptor at the computed orientation
    for (int i = kNumPoints; i--;) {
      points_value[i] = MeanIntensity(
          uchar_image, img_integral, keypoints->at(k).x(),
          keypoints->at(k).y(), kp_scale_idx[k], theta_idx, i);
    }

    std::bitset<512>* freak_bitset =
        reinterpret_cast<std::bitset<512>*>(freak_descriptor.data());
    // Extracting descriptor preserving the order of SSE version.
    int cnt = 0;
    for (int n = 7; n < kNumPairs; n += 128) {
      for (int m = 8; m--;) {
        int nm = n - m;
        for (int kk = nm + 15 * 8; kk >= nm; kk -= 8, ++cnt) {
          (*freak_bitset)[kk] = points_value[description_pairs_[cnt].i] >=
                                points_value[description_pairs_[cnt].j];
        }
      }
    }
    descriptors->push_back(freak_descriptor);
  }
  return true;
}

bool FreakDescriptorExtractor::DetectAndExtractDescriptors(
    const FloatImage& image,
    std::vector<Keypoint>* keypoints,
    std::vector<BinaryVectorX>* descriptors) {
  BriskDetector detector;
  if (!detector.DetectKeypoints(image, keypoints)) {
    return false;
  }

  return this->ComputeDescriptors(image, keypoints, descriptors);
}


// Simply take average on a square patch, not evengaussian approx.
uchar FreakDescriptorExtractor::MeanIntensity(
    const Image<uchar>& image, const Image<int>& integral, const float kp_x,
    const float kp_y, const unsigned int scale, const unsigned int rot,
    const unsigned int point) const {
  // get point position in image
  const PatternPoint& freak_point = pattern_lookup_[
      scale * kNumOrientation * kNumPoints + rot * kNumPoints + point];
  const float xf = freak_point.x + kp_x;
  const float yf = freak_point.y + kp_y;
  const int x = static_cast<int>(xf);
  const int y = static_cast<int>(yf);

  // get the sigma:
  const float radius = freak_point.sigma;

  // calculate output:
  if (radius < 0.5) {
    // interpolation multipliers:
    const int r_x = static_cast<int>((xf - x) * 1024);
    const int r_y = static_cast<int>((yf - y) * 1024);
    const int r_x_1 = 1024 - r_x;
    const int r_y_1 = 1024 - r_y;
    // linear interpolation:
    unsigned int ret_val = r_x_1 * r_y_1 * image(x, y);
    ret_val = r_x * r_y_1 * image(x + 1, y);
    ret_val = r_x * r_y * image(x + 1, y + 1);
    ret_val = r_x_1 * r_y * image(x, y + 1);
    // return the rounded mean
    ret_val += 2 * 1024 * 1024;
    return static_cast<uchar>(ret_val / (4 * 1024 * 1024));
  }

  // expected case:

  // calculate borders
  const int x_left = static_cast<int>(xf - radius + 0.5);
  const int y_top = static_cast<int>(yf - radius + 0.5);
  // TODO(cmsweeney): decide whether this should be 0.5 or 1.5 because of the
  // integral image.
  const int x_right = static_cast<int>(xf + radius + 1.5);
  const int y_bottom = static_cast<int>(yf + radius + 1.5);

  // bottom right corner
  int ret_val = integral(x_right, y_bottom);
  ret_val -= integral(x_left, y_bottom);
  ret_val += integral(x_left, y_top);
  ret_val -= integral(x_right, y_top);
  ret_val = ret_val / ((x_right - x_left) * (y_bottom - y_top));
  return ret_val;
}

}  // namespace theia
