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

#include "theia/sfm/pose/four_point_focal_length_helper.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <glog/logging.h>
#include <vector>

namespace theia {

using Eigen::Matrix;
using Eigen::Vector3d;

// Helper method to the P4Pf algorithm that computes the grobner basis that
// allows for a solution to the problem.
bool FourPointFocalLengthHelper(
    const double glab, const double glac, const double glad, const double glbc,
    const double glbd, const double glcd,
    const Eigen::Matrix<double, 2, 4>& features_normalized,
    std::vector<double>* f, std::vector<Eigen::Vector3d>* depths) {
  // Set convenience vars.
  const double a1 = features_normalized(0, 0);
  const double a2 = features_normalized(1, 0);
  const double b1 = features_normalized(0, 1);
  const double b2 = features_normalized(1, 1);
  const double c1 = features_normalized(0, 2);
  const double c2 = features_normalized(1, 2);
  const double d1 = features_normalized(0, 3);
  const double d2 = features_normalized(1, 3);

  // Precompute coefficients.
  Matrix<double, 36, 1> coeffs;
  coeffs(0) = 1.0;
  coeffs(1) = (glab * (-1.0 / 2.0)) / glad - (glac * (1.0 / 2.0)) / glad +
              (glbc * (1.0 / 2.0)) / glad;
  coeffs(2) = -1.0;
  coeffs(3) = b1 * c1 + b2 * c2;
  coeffs(4) = glab / glad + glac / glad - glbc / glad;
  coeffs(5) = ((d1 * d1) * glab * (-1.0 / 2.0)) / glad -
              ((d1 * d1) * glac * (1.0 / 2.0)) / glad -
              ((d2 * d2) * glab * (1.0 / 2.0)) / glad -
              ((d2 * d2) * glac * (1.0 / 2.0)) / glad +
              ((d1 * d1) * glbc * (1.0 / 2.0)) / glad +
              ((d2 * d2) * glbc * (1.0 / 2.0)) / glad;
  coeffs(6) = (glab * (-1.0 / 2.0)) / glad - (glac * (1.0 / 2.0)) / glad +
              (glbc * (1.0 / 2.0)) / glad + 1.0;
  coeffs(7) = -a1 * b1 - a2 * b2;
  coeffs(8) = -a1 * c1 - a2 * c2;
  coeffs(9) = (a1 * d1 * glab) / glad + (a1 * d1 * glac) / glad +
              (a2 * d2 * glab) / glad + (a2 * d2 * glac) / glad -
              (a1 * d1 * glbc) / glad - (a2 * d2 * glbc) / glad;
  coeffs(10) = a1 * a1 + a2 * a2 - ((a1 * a1) * glab * (1.0 / 2.0)) / glad -
               ((a1 * a1) * glac * (1.0 / 2.0)) / glad -
               ((a2 * a2) * glab * (1.0 / 2.0)) / glad -
               ((a2 * a2) * glac * (1.0 / 2.0)) / glad +
               ((a1 * a1) * glbc * (1.0 / 2.0)) / glad +
               ((a2 * a2) * glbc * (1.0 / 2.0)) / glad;
  coeffs(11) = -glac / glad;
  coeffs(12) = -2.0;
  coeffs(13) = c1 * c1 + c2 * c2;
  coeffs(14) = (glac * 2.0) / glad;
  coeffs(15) = -((d1 * d1) * glac) / glad - ((d2 * d2) * glac) / glad;
  coeffs(16) = -glac / glad + 1.0;
  coeffs(17) = a1 * c1 * -2.0 - a2 * c2 * 2.0;
  coeffs(18) = (a1 * d1 * glac * 2.0) / glad + (a2 * d2 * glac * 2.0) / glad;
  coeffs(19) =
      a1 * a1 + a2 * a2 - ((a1 * a1) * glac) / glad - ((a2 * a2) * glac) / glad;
  coeffs(20) =
      (glab * (-1.0 / 2.0)) / glad + (glbd * (1.0 / 2.0)) / glad - 1.0 / 2.0;
  coeffs(21) = glab / glad - glbd / glad;
  coeffs(22) = b1 * d1 + b2 * d2;
  coeffs(23) = (d1 * d1) * (-1.0 / 2.0) - (d2 * d2) * (1.0 / 2.0) -
               ((d1 * d1) * glab * (1.0 / 2.0)) / glad -
               ((d2 * d2) * glab * (1.0 / 2.0)) / glad +
               ((d1 * d1) * glbd * (1.0 / 2.0)) / glad +
               ((d2 * d2) * glbd * (1.0 / 2.0)) / glad;
  coeffs(24) =
      (glab * (-1.0 / 2.0)) / glad + (glbd * (1.0 / 2.0)) / glad + 1.0 / 2.0;
  coeffs(25) = -a1 * b1 - a2 * b2;
  coeffs(26) = (a1 * d1 * glab) / glad + (a2 * d2 * glab) / glad -
               (a1 * d1 * glbd) / glad - (a2 * d2 * glbd) / glad;
  coeffs(27) = (a1 * a1) * (1.0 / 2.0) + (a2 * a2) * (1.0 / 2.0) -
               ((a1 * a1) * glab * (1.0 / 2.0)) / glad -
               ((a2 * a2) * glab * (1.0 / 2.0)) / glad +
               ((a1 * a1) * glbd * (1.0 / 2.0)) / glad +
               ((a2 * a2) * glbd * (1.0 / 2.0)) / glad;
  coeffs(28) =
      (glac * (-1.0 / 2.0)) / glad + (glcd * (1.0 / 2.0)) / glad - 1.0 / 2.0;
  coeffs(29) = glac / glad - glcd / glad;
  coeffs(30) = c1 * d1 + c2 * d2;
  coeffs(31) = (d1 * d1) * (-1.0 / 2.0) - (d2 * d2) * (1.0 / 2.0) -
               ((d1 * d1) * glac * (1.0 / 2.0)) / glad -
               ((d2 * d2) * glac * (1.0 / 2.0)) / glad +
               ((d1 * d1) * glcd * (1.0 / 2.0)) / glad +
               ((d2 * d2) * glcd * (1.0 / 2.0)) / glad;
  coeffs(32) =
      (glac * (-1.0 / 2.0)) / glad + (glcd * (1.0 / 2.0)) / glad + 1.0 / 2.0;
  coeffs(33) = -a1 * c1 - a2 * c2;
  coeffs(34) = (a1 * d1 * glac) / glad + (a2 * d2 * glac) / glad -
               (a1 * d1 * glcd) / glad - (a2 * d2 * glcd) / glad;
  coeffs(35) = (a1 * a1) * (1.0 / 2.0) + (a2 * a2) * (1.0 / 2.0) -
               ((a1 * a1) * glac * (1.0 / 2.0)) / glad -
               ((a2 * a2) * glac * (1.0 / 2.0)) / glad +
               ((a1 * a1) * glcd * (1.0 / 2.0)) / glad +
               ((a2 * a2) * glcd * (1.0 / 2.0)) / glad;

  Matrix<double, 78, 78> coefficient_matrix;
  coefficient_matrix.setZero();
  Matrix<double, 78, 10> b;
  b.setZero();

  // Set matrix coefficients. Slightly awkward because of indexing differences
  // from matlab.
  coefficient_matrix(0, 64) = coeffs(1 - 1);
  coefficient_matrix(0, 71) = coeffs(2 - 1);
  coefficient_matrix(0, 75) = coeffs(3 - 1);
  coefficient_matrix(0, 76) = coeffs(3 - 1);
  coefficient_matrix(0, 77) = coeffs(4 - 1);
  coefficient_matrix(1, 65) = coeffs(1 - 1);
  coefficient_matrix(1, 71) = coeffs(12 - 1);
  coefficient_matrix(1, 76) = coeffs(13 - 1);
  coefficient_matrix(2, 67) = coeffs(1 - 1);
  coefficient_matrix(2, 71) = coeffs(21 - 1);
  coefficient_matrix(2, 75) = coeffs(3 - 1);
  coefficient_matrix(3, 68) = coeffs(1 - 1);
  coefficient_matrix(3, 71) = coeffs(29 - 1);
  coefficient_matrix(3, 76) = coeffs(3 - 1);
  coefficient_matrix(4, 51) = coeffs(1 - 1);
  coefficient_matrix(4, 59) = coeffs(2 - 1);
  coefficient_matrix(4, 67) = coeffs(3 - 1);
  coefficient_matrix(4, 68) = coeffs(3 - 1);
  coefficient_matrix(4, 69) = coeffs(4 - 1);
  coefficient_matrix(4, 71) = coeffs(5 - 1);
  coefficient_matrix(4, 74) = coeffs(6 - 1);
  coefficient_matrix(5, 46) = coeffs(1 - 1);
  coefficient_matrix(5, 56) = coeffs(2 - 1);
  coefficient_matrix(5, 64) = coeffs(3 - 1);
  coefficient_matrix(5, 65) = coeffs(3 - 1);
  coefficient_matrix(5, 66) = coeffs(4 - 1);
  coefficient_matrix(5, 68) = coeffs(5 - 1);
  coefficient_matrix(5, 73) = coeffs(6 - 1);
  coefficient_matrix(5, 76) = coeffs(7 - 1);
  coefficient_matrix(5, 77) = coeffs(8 - 1);
  coefficient_matrix(6, 44) = coeffs(1 - 1);
  coefficient_matrix(6, 54) = coeffs(2 - 1);
  coefficient_matrix(6, 62) = coeffs(3 - 1);
  coefficient_matrix(6, 63) = coeffs(3 - 1);
  coefficient_matrix(6, 64) = coeffs(4 - 1);
  coefficient_matrix(6, 71) = coeffs(6 - 1);
  coefficient_matrix(6, 75) = coeffs(8 - 1);
  coefficient_matrix(6, 76) = coeffs(9 - 1);
  coefficient_matrix(7, 52) = coeffs(1 - 1);
  coefficient_matrix(7, 59) = coeffs(12 - 1);
  coefficient_matrix(7, 68) = coeffs(13 - 1);
  coefficient_matrix(7, 70) = coeffs(14 - 1);
  coefficient_matrix(7, 71) = coeffs(15 - 1);
  coefficient_matrix(7, 74) = coeffs(16 - 1);
  coefficient_matrix(8, 46) = coeffs(1 - 1);
  coefficient_matrix(8, 55) = coeffs(12 - 1);
  coefficient_matrix(8, 64) = coeffs(13 - 1);
  coefficient_matrix(8, 66) = coeffs(14 - 1);
  coefficient_matrix(8, 67) = coeffs(15 - 1);
  coefficient_matrix(8, 72) = coeffs(16 - 1);
  coefficient_matrix(8, 75) = coeffs(17 - 1);
  coefficient_matrix(8, 77) = coeffs(18 - 1);
  coefficient_matrix(9, 45) = coeffs(1 - 1);
  coefficient_matrix(9, 54) = coeffs(12 - 1);
  coefficient_matrix(9, 63) = coeffs(13 - 1);
  coefficient_matrix(9, 65) = coeffs(14 - 1);
  coefficient_matrix(9, 71) = coeffs(16 - 1);
  coefficient_matrix(9, 76) = coeffs(18 - 1);
  coefficient_matrix(10, 55) = coeffs(1 - 1);
  coefficient_matrix(10, 59) = coeffs(21 - 1);
  coefficient_matrix(10, 67) = coeffs(3 - 1);
  coefficient_matrix(10, 71) = coeffs(22 - 1);
  coefficient_matrix(10, 72) = coeffs(23 - 1);
  coefficient_matrix(10, 74) = coeffs(24 - 1);
  coefficient_matrix(11, 51) = coeffs(1 - 1);
  coefficient_matrix(11, 56) = coeffs(21 - 1);
  coefficient_matrix(11, 64) = coeffs(3 - 1);
  coefficient_matrix(11, 68) = coeffs(22 - 1);
  coefficient_matrix(11, 69) = coeffs(23 - 1);
  coefficient_matrix(11, 73) = coeffs(24 - 1);
  coefficient_matrix(11, 76) = coeffs(25 - 1);
  coefficient_matrix(11, 77) = coeffs(26 - 1);
  coefficient_matrix(12, 49) = coeffs(1 - 1);
  coefficient_matrix(12, 54) = coeffs(21 - 1);
  coefficient_matrix(12, 62) = coeffs(3 - 1);
  coefficient_matrix(12, 67) = coeffs(23 - 1);
  coefficient_matrix(12, 71) = coeffs(24 - 1);
  coefficient_matrix(12, 75) = coeffs(26 - 1);
  coefficient_matrix(13, 56) = coeffs(1 - 1);
  coefficient_matrix(13, 59) = coeffs(29 - 1);
  coefficient_matrix(13, 68) = coeffs(3 - 1);
  coefficient_matrix(13, 71) = coeffs(30 - 1);
  coefficient_matrix(13, 73) = coeffs(31 - 1);
  coefficient_matrix(13, 74) = coeffs(32 - 1);
  coefficient_matrix(14, 52) = coeffs(1 - 1);
  coefficient_matrix(14, 56) = coeffs(29 - 1);
  coefficient_matrix(14, 65) = coeffs(3 - 1);
  coefficient_matrix(14, 68) = coeffs(30 - 1);
  coefficient_matrix(14, 70) = coeffs(31 - 1);
  coefficient_matrix(14, 73) = coeffs(32 - 1);
  coefficient_matrix(14, 76) = coeffs(33 - 1);
  coefficient_matrix(15, 51) = coeffs(1 - 1);
  coefficient_matrix(15, 55) = coeffs(29 - 1);
  coefficient_matrix(15, 64) = coeffs(3 - 1);
  coefficient_matrix(15, 67) = coeffs(30 - 1);
  coefficient_matrix(15, 69) = coeffs(31 - 1);
  coefficient_matrix(15, 72) = coeffs(32 - 1);
  coefficient_matrix(15, 75) = coeffs(33 - 1);
  coefficient_matrix(15, 77) = coeffs(34 - 1);
  coefficient_matrix(16, 50) = coeffs(1 - 1);
  coefficient_matrix(16, 54) = coeffs(29 - 1);
  coefficient_matrix(16, 63) = coeffs(3 - 1);
  coefficient_matrix(16, 68) = coeffs(31 - 1);
  coefficient_matrix(16, 71) = coeffs(32 - 1);
  coefficient_matrix(16, 76) = coeffs(34 - 1);
  coefficient_matrix(17, 37) = coeffs(1 - 1);
  coefficient_matrix(17, 43) = coeffs(2 - 1);
  coefficient_matrix(17, 55) = coeffs(3 - 1);
  coefficient_matrix(17, 56) = coeffs(3 - 1);
  coefficient_matrix(17, 57) = coeffs(4 - 1);
  coefficient_matrix(17, 59) = coeffs(5 - 1);
  coefficient_matrix(17, 71) = coeffs(7 - 1);
  coefficient_matrix(17, 72) = coeffs(8 - 1);
  coefficient_matrix(17, 73) = coeffs(9 - 1);
  coefficient_matrix(17, 74) = coeffs(10 - 1);
  coefficient_matrix(18, 32) = coeffs(1 - 1);
  coefficient_matrix(18, 42) = coeffs(2 - 1);
  coefficient_matrix(18, 51) = coeffs(3 - 1);
  coefficient_matrix(18, 52) = coeffs(3 - 1);
  coefficient_matrix(18, 53) = coeffs(4 - 1);
  coefficient_matrix(18, 56) = coeffs(5 - 1);
  coefficient_matrix(18, 61) = coeffs(6 - 1);
  coefficient_matrix(18, 68) = coeffs(7 - 1);
  coefficient_matrix(18, 69) = coeffs(8 - 1);
  coefficient_matrix(18, 70) = coeffs(9 - 1);
  coefficient_matrix(18, 73) = coeffs(10 - 1);
  coefficient_matrix(19, 30) = coeffs(1 - 1);
  coefficient_matrix(19, 40) = coeffs(2 - 1);
  coefficient_matrix(19, 49) = coeffs(3 - 1);
  coefficient_matrix(19, 50) = coeffs(3 - 1);
  coefficient_matrix(19, 51) = coeffs(4 - 1);
  coefficient_matrix(19, 54) = coeffs(5 - 1);
  coefficient_matrix(19, 59) = coeffs(6 - 1);
  coefficient_matrix(19, 67) = coeffs(8 - 1);
  coefficient_matrix(19, 68) = coeffs(9 - 1);
  coefficient_matrix(19, 71) = coeffs(10 - 1);
  coefficient_matrix(20, 27) = coeffs(1 - 1);
  coefficient_matrix(20, 38) = coeffs(2 - 1);
  coefficient_matrix(20, 46) = coeffs(3 - 1);
  coefficient_matrix(20, 47) = coeffs(3 - 1);
  coefficient_matrix(20, 48) = coeffs(4 - 1);
  coefficient_matrix(20, 52) = coeffs(5 - 1);
  coefficient_matrix(20, 58) = coeffs(6 - 1);
  coefficient_matrix(20, 65) = coeffs(7 - 1);
  coefficient_matrix(20, 66) = coeffs(8 - 1);
  coefficient_matrix(20, 70) = coeffs(10 - 1);
  coefficient_matrix(21, 26) = coeffs(1 - 1);
  coefficient_matrix(21, 36) = coeffs(2 - 1);
  coefficient_matrix(21, 44) = coeffs(3 - 1);
  coefficient_matrix(21, 45) = coeffs(3 - 1);
  coefficient_matrix(21, 46) = coeffs(4 - 1);
  coefficient_matrix(21, 50) = coeffs(5 - 1);
  coefficient_matrix(21, 56) = coeffs(6 - 1);
  coefficient_matrix(21, 63) = coeffs(7 - 1);
  coefficient_matrix(21, 64) = coeffs(8 - 1);
  coefficient_matrix(21, 65) = coeffs(9 - 1);
  coefficient_matrix(21, 68) = coeffs(10 - 1);
  coefficient_matrix(21, 76) = coeffs(11 - 1);
  coefficient_matrix(22, 24) = coeffs(1 - 1);
  coefficient_matrix(22, 44) = coeffs(4 - 1);
  coefficient_matrix(22, 54) = coeffs(6 - 1);
  coefficient_matrix(22, 62) = coeffs(8 - 1);
  coefficient_matrix(22, 63) = coeffs(9 - 1);
  coefficient_matrix(23, 38) = coeffs(1 - 1);
  coefficient_matrix(23, 43) = coeffs(12 - 1);
  coefficient_matrix(23, 56) = coeffs(13 - 1);
  coefficient_matrix(23, 58) = coeffs(14 - 1);
  coefficient_matrix(23, 59) = coeffs(15 - 1);
  coefficient_matrix(23, 71) = coeffs(17 - 1);
  coefficient_matrix(23, 73) = coeffs(18 - 1);
  coefficient_matrix(23, 74) = coeffs(19 - 1);
  coefficient_matrix(24, 33) = coeffs(1 - 1);
  coefficient_matrix(24, 42) = coeffs(12 - 1);
  coefficient_matrix(24, 52) = coeffs(13 - 1);
  coefficient_matrix(24, 56) = coeffs(15 - 1);
  coefficient_matrix(24, 61) = coeffs(16 - 1);
  coefficient_matrix(24, 68) = coeffs(17 - 1);
  coefficient_matrix(24, 70) = coeffs(18 - 1);
  coefficient_matrix(24, 73) = coeffs(19 - 1);
  coefficient_matrix(25, 32) = coeffs(1 - 1);
  coefficient_matrix(25, 41) = coeffs(12 - 1);
  coefficient_matrix(25, 51) = coeffs(13 - 1);
  coefficient_matrix(25, 53) = coeffs(14 - 1);
  coefficient_matrix(25, 55) = coeffs(15 - 1);
  coefficient_matrix(25, 60) = coeffs(16 - 1);
  coefficient_matrix(25, 67) = coeffs(17 - 1);
  coefficient_matrix(25, 69) = coeffs(18 - 1);
  coefficient_matrix(25, 72) = coeffs(19 - 1);
  coefficient_matrix(26, 31) = coeffs(1 - 1);
  coefficient_matrix(26, 40) = coeffs(12 - 1);
  coefficient_matrix(26, 50) = coeffs(13 - 1);
  coefficient_matrix(26, 52) = coeffs(14 - 1);
  coefficient_matrix(26, 54) = coeffs(15 - 1);
  coefficient_matrix(26, 59) = coeffs(16 - 1);
  coefficient_matrix(26, 68) = coeffs(18 - 1);
  coefficient_matrix(26, 71) = coeffs(19 - 1);
  coefficient_matrix(27, 27) = coeffs(1 - 1);
  coefficient_matrix(27, 37) = coeffs(12 - 1);
  coefficient_matrix(27, 46) = coeffs(13 - 1);
  coefficient_matrix(27, 48) = coeffs(14 - 1);
  coefficient_matrix(27, 51) = coeffs(15 - 1);
  coefficient_matrix(27, 57) = coeffs(16 - 1);
  coefficient_matrix(27, 64) = coeffs(17 - 1);
  coefficient_matrix(27, 66) = coeffs(18 - 1);
  coefficient_matrix(27, 69) = coeffs(19 - 1);
  coefficient_matrix(27, 77) = coeffs(20 - 1);
  coefficient_matrix(28, 26) = coeffs(1 - 1);
  coefficient_matrix(28, 35) = coeffs(12 - 1);
  coefficient_matrix(28, 44) = coeffs(13 - 1);
  coefficient_matrix(28, 46) = coeffs(14 - 1);
  coefficient_matrix(28, 49) = coeffs(15 - 1);
  coefficient_matrix(28, 55) = coeffs(16 - 1);
  coefficient_matrix(28, 62) = coeffs(17 - 1);
  coefficient_matrix(28, 64) = coeffs(18 - 1);
  coefficient_matrix(28, 67) = coeffs(19 - 1);
  coefficient_matrix(28, 75) = coeffs(20 - 1);
  coefficient_matrix(29, 25) = coeffs(1 - 1);
  coefficient_matrix(29, 45) = coeffs(14 - 1);
  coefficient_matrix(29, 54) = coeffs(16 - 1);
  coefficient_matrix(29, 63) = coeffs(18 - 1);
  coefficient_matrix(30, 41) = coeffs(1 - 1);
  coefficient_matrix(30, 43) = coeffs(21 - 1);
  coefficient_matrix(30, 55) = coeffs(3 - 1);
  coefficient_matrix(30, 59) = coeffs(22 - 1);
  coefficient_matrix(30, 60) = coeffs(23 - 1);
  coefficient_matrix(30, 71) = coeffs(25 - 1);
  coefficient_matrix(30, 72) = coeffs(26 - 1);
  coefficient_matrix(30, 74) = coeffs(27 - 1);
  coefficient_matrix(31, 37) = coeffs(1 - 1);
  coefficient_matrix(31, 42) = coeffs(21 - 1);
  coefficient_matrix(31, 51) = coeffs(3 - 1);
  coefficient_matrix(31, 56) = coeffs(22 - 1);
  coefficient_matrix(31, 57) = coeffs(23 - 1);
  coefficient_matrix(31, 61) = coeffs(24 - 1);
  coefficient_matrix(31, 68) = coeffs(25 - 1);
  coefficient_matrix(31, 69) = coeffs(26 - 1);
  coefficient_matrix(31, 73) = coeffs(27 - 1);
  coefficient_matrix(32, 35) = coeffs(1 - 1);
  coefficient_matrix(32, 40) = coeffs(21 - 1);
  coefficient_matrix(32, 49) = coeffs(3 - 1);
  coefficient_matrix(32, 54) = coeffs(22 - 1);
  coefficient_matrix(32, 55) = coeffs(23 - 1);
  coefficient_matrix(32, 59) = coeffs(24 - 1);
  coefficient_matrix(32, 67) = coeffs(26 - 1);
  coefficient_matrix(32, 71) = coeffs(27 - 1);
  coefficient_matrix(33, 32) = coeffs(1 - 1);
  coefficient_matrix(33, 38) = coeffs(21 - 1);
  coefficient_matrix(33, 46) = coeffs(3 - 1);
  coefficient_matrix(33, 52) = coeffs(22 - 1);
  coefficient_matrix(33, 53) = coeffs(23 - 1);
  coefficient_matrix(33, 58) = coeffs(24 - 1);
  coefficient_matrix(33, 65) = coeffs(25 - 1);
  coefficient_matrix(33, 66) = coeffs(26 - 1);
  coefficient_matrix(33, 70) = coeffs(27 - 1);
  coefficient_matrix(34, 30) = coeffs(1 - 1);
  coefficient_matrix(34, 36) = coeffs(21 - 1);
  coefficient_matrix(34, 44) = coeffs(3 - 1);
  coefficient_matrix(34, 50) = coeffs(22 - 1);
  coefficient_matrix(34, 51) = coeffs(23 - 1);
  coefficient_matrix(34, 56) = coeffs(24 - 1);
  coefficient_matrix(34, 63) = coeffs(25 - 1);
  coefficient_matrix(34, 64) = coeffs(26 - 1);
  coefficient_matrix(34, 68) = coeffs(27 - 1);
  coefficient_matrix(34, 76) = coeffs(28 - 1);
  coefficient_matrix(35, 28) = coeffs(1 - 1);
  coefficient_matrix(35, 49) = coeffs(23 - 1);
  coefficient_matrix(35, 54) = coeffs(24 - 1);
  coefficient_matrix(35, 62) = coeffs(26 - 1);
  coefficient_matrix(36, 42) = coeffs(1 - 1);
  coefficient_matrix(36, 43) = coeffs(29 - 1);
  coefficient_matrix(36, 56) = coeffs(3 - 1);
  coefficient_matrix(36, 59) = coeffs(30 - 1);
  coefficient_matrix(36, 61) = coeffs(31 - 1);
  coefficient_matrix(36, 71) = coeffs(33 - 1);
  coefficient_matrix(36, 73) = coeffs(34 - 1);
  coefficient_matrix(36, 74) = coeffs(35 - 1);
  coefficient_matrix(37, 38) = coeffs(1 - 1);
  coefficient_matrix(37, 42) = coeffs(29 - 1);
  coefficient_matrix(37, 52) = coeffs(3 - 1);
  coefficient_matrix(37, 56) = coeffs(30 - 1);
  coefficient_matrix(37, 58) = coeffs(31 - 1);
  coefficient_matrix(37, 61) = coeffs(32 - 1);
  coefficient_matrix(37, 68) = coeffs(33 - 1);
  coefficient_matrix(37, 70) = coeffs(34 - 1);
  coefficient_matrix(37, 73) = coeffs(35 - 1);
  coefficient_matrix(38, 37) = coeffs(1 - 1);
  coefficient_matrix(38, 41) = coeffs(29 - 1);
  coefficient_matrix(38, 51) = coeffs(3 - 1);
  coefficient_matrix(38, 55) = coeffs(30 - 1);
  coefficient_matrix(38, 57) = coeffs(31 - 1);
  coefficient_matrix(38, 60) = coeffs(32 - 1);
  coefficient_matrix(38, 67) = coeffs(33 - 1);
  coefficient_matrix(38, 69) = coeffs(34 - 1);
  coefficient_matrix(38, 72) = coeffs(35 - 1);
  coefficient_matrix(39, 36) = coeffs(1 - 1);
  coefficient_matrix(39, 40) = coeffs(29 - 1);
  coefficient_matrix(39, 50) = coeffs(3 - 1);
  coefficient_matrix(39, 54) = coeffs(30 - 1);
  coefficient_matrix(39, 56) = coeffs(31 - 1);
  coefficient_matrix(39, 59) = coeffs(32 - 1);
  coefficient_matrix(39, 68) = coeffs(34 - 1);
  coefficient_matrix(39, 71) = coeffs(35 - 1);
  coefficient_matrix(40, 33) = coeffs(1 - 1);
  coefficient_matrix(40, 38) = coeffs(29 - 1);
  coefficient_matrix(40, 47) = coeffs(3 - 1);
  coefficient_matrix(40, 52) = coeffs(30 - 1);
  coefficient_matrix(40, 58) = coeffs(32 - 1);
  coefficient_matrix(40, 65) = coeffs(33 - 1);
  coefficient_matrix(40, 70) = coeffs(35 - 1);
  coefficient_matrix(41, 32) = coeffs(1 - 1);
  coefficient_matrix(41, 37) = coeffs(29 - 1);
  coefficient_matrix(41, 46) = coeffs(3 - 1);
  coefficient_matrix(41, 51) = coeffs(30 - 1);
  coefficient_matrix(41, 53) = coeffs(31 - 1);
  coefficient_matrix(41, 57) = coeffs(32 - 1);
  coefficient_matrix(41, 64) = coeffs(33 - 1);
  coefficient_matrix(41, 66) = coeffs(34 - 1);
  coefficient_matrix(41, 69) = coeffs(35 - 1);
  coefficient_matrix(41, 77) = coeffs(36 - 1);
  coefficient_matrix(42, 31) = coeffs(1 - 1);
  coefficient_matrix(42, 36) = coeffs(29 - 1);
  coefficient_matrix(42, 45) = coeffs(3 - 1);
  coefficient_matrix(42, 50) = coeffs(30 - 1);
  coefficient_matrix(42, 52) = coeffs(31 - 1);
  coefficient_matrix(42, 56) = coeffs(32 - 1);
  coefficient_matrix(42, 63) = coeffs(33 - 1);
  coefficient_matrix(42, 65) = coeffs(34 - 1);
  coefficient_matrix(42, 68) = coeffs(35 - 1);
  coefficient_matrix(42, 76) = coeffs(36 - 1);
  coefficient_matrix(43, 30) = coeffs(1 - 1);
  coefficient_matrix(43, 35) = coeffs(29 - 1);
  coefficient_matrix(43, 44) = coeffs(3 - 1);
  coefficient_matrix(43, 49) = coeffs(30 - 1);
  coefficient_matrix(43, 51) = coeffs(31 - 1);
  coefficient_matrix(43, 55) = coeffs(32 - 1);
  coefficient_matrix(43, 62) = coeffs(33 - 1);
  coefficient_matrix(43, 64) = coeffs(34 - 1);
  coefficient_matrix(43, 67) = coeffs(35 - 1);
  coefficient_matrix(43, 75) = coeffs(36 - 1);
  coefficient_matrix(44, 29) = coeffs(1 - 1);
  coefficient_matrix(44, 50) = coeffs(31 - 1);
  coefficient_matrix(44, 54) = coeffs(32 - 1);
  coefficient_matrix(44, 63) = coeffs(34 - 1);
  coefficient_matrix(45, 19) = coeffs(1 - 1);
  coefficient_matrix(45, 41) = coeffs(3 - 1);
  coefficient_matrix(45, 42) = coeffs(3 - 1);
  coefficient_matrix(45, 43) = coeffs(5 - 1);
  coefficient_matrix(45, 59) = coeffs(7 - 1);
  coefficient_matrix(45, 60) = coeffs(8 - 1);
  coefficient_matrix(45, 61) = coeffs(9 - 1);
  coefficient_matrix(45, 74) = coeffs(11 - 1);
  coefficient_matrix(46, 15) = coeffs(1 - 1);
  coefficient_matrix(46, 23) = coeffs(2 - 1);
  coefficient_matrix(46, 37) = coeffs(3 - 1);
  coefficient_matrix(46, 38) = coeffs(3 - 1);
  coefficient_matrix(46, 39) = coeffs(4 - 1);
  coefficient_matrix(46, 42) = coeffs(5 - 1);
  coefficient_matrix(46, 56) = coeffs(7 - 1);
  coefficient_matrix(46, 57) = coeffs(8 - 1);
  coefficient_matrix(46, 58) = coeffs(9 - 1);
  coefficient_matrix(46, 61) = coeffs(10 - 1);
  coefficient_matrix(46, 73) = coeffs(11 - 1);
  coefficient_matrix(47, 13) = coeffs(1 - 1);
  coefficient_matrix(47, 21) = coeffs(2 - 1);
  coefficient_matrix(47, 35) = coeffs(3 - 1);
  coefficient_matrix(47, 36) = coeffs(3 - 1);
  coefficient_matrix(47, 37) = coeffs(4 - 1);
  coefficient_matrix(47, 40) = coeffs(5 - 1);
  coefficient_matrix(47, 43) = coeffs(6 - 1);
  coefficient_matrix(47, 54) = coeffs(7 - 1);
  coefficient_matrix(47, 55) = coeffs(8 - 1);
  coefficient_matrix(47, 56) = coeffs(9 - 1);
  coefficient_matrix(47, 59) = coeffs(10 - 1);
  coefficient_matrix(47, 71) = coeffs(11 - 1);
  coefficient_matrix(48, 10) = coeffs(1 - 1);
  coefficient_matrix(48, 20) = coeffs(2 - 1);
  coefficient_matrix(48, 32) = coeffs(3 - 1);
  coefficient_matrix(48, 33) = coeffs(3 - 1);
  coefficient_matrix(48, 34) = coeffs(4 - 1);
  coefficient_matrix(48, 38) = coeffs(5 - 1);
  coefficient_matrix(48, 52) = coeffs(7 - 1);
  coefficient_matrix(48, 53) = coeffs(8 - 1);
  coefficient_matrix(48, 58) = coeffs(10 - 1);
  coefficient_matrix(48, 70) = coeffs(11 - 1);
  coefficient_matrix(49, 9) = coeffs(1 - 1);
  coefficient_matrix(49, 18) = coeffs(2 - 1);
  coefficient_matrix(49, 30) = coeffs(3 - 1);
  coefficient_matrix(49, 31) = coeffs(3 - 1);
  coefficient_matrix(49, 32) = coeffs(4 - 1);
  coefficient_matrix(49, 36) = coeffs(5 - 1);
  coefficient_matrix(49, 42) = coeffs(6 - 1);
  coefficient_matrix(49, 50) = coeffs(7 - 1);
  coefficient_matrix(49, 51) = coeffs(8 - 1);
  coefficient_matrix(49, 52) = coeffs(9 - 1);
  coefficient_matrix(49, 56) = coeffs(10 - 1);
  coefficient_matrix(49, 68) = coeffs(11 - 1);
  coefficient_matrix(50, 7) = coeffs(1 - 1);
  coefficient_matrix(50, 16) = coeffs(2 - 1);
  coefficient_matrix(50, 28) = coeffs(3 - 1);
  coefficient_matrix(50, 29) = coeffs(3 - 1);
  coefficient_matrix(50, 30) = coeffs(4 - 1);
  coefficient_matrix(50, 40) = coeffs(6 - 1);
  coefficient_matrix(50, 49) = coeffs(8 - 1);
  coefficient_matrix(50, 50) = coeffs(9 - 1);
  coefficient_matrix(50, 54) = coeffs(10 - 1);
  coefficient_matrix(51, 6) = coeffs(1 - 1);
  coefficient_matrix(51, 12) = coeffs(2 - 1);
  coefficient_matrix(51, 24) = coeffs(3 - 1);
  coefficient_matrix(51, 25) = coeffs(3 - 1);
  coefficient_matrix(51, 26) = coeffs(4 - 1);
  coefficient_matrix(51, 29) = coeffs(5 - 1);
  coefficient_matrix(51, 36) = coeffs(6 - 1);
  coefficient_matrix(51, 44) = coeffs(8 - 1);
  coefficient_matrix(51, 45) = coeffs(9 - 1);
  coefficient_matrix(51, 50) = coeffs(10 - 1);
  coefficient_matrix(51, 63) = coeffs(11 - 1);
  coefficient_matrix(52, 20) = coeffs(1 - 1);
  coefficient_matrix(52, 42) = coeffs(13 - 1);
  coefficient_matrix(52, 43) = coeffs(15 - 1);
  coefficient_matrix(52, 59) = coeffs(17 - 1);
  coefficient_matrix(52, 61) = coeffs(18 - 1);
  coefficient_matrix(52, 74) = coeffs(20 - 1);
  coefficient_matrix(53, 15) = coeffs(1 - 1);
  coefficient_matrix(53, 22) = coeffs(12 - 1);
  coefficient_matrix(53, 37) = coeffs(13 - 1);
  coefficient_matrix(53, 39) = coeffs(14 - 1);
  coefficient_matrix(53, 41) = coeffs(15 - 1);
  coefficient_matrix(53, 55) = coeffs(17 - 1);
  coefficient_matrix(53, 57) = coeffs(18 - 1);
  coefficient_matrix(53, 60) = coeffs(19 - 1);
  coefficient_matrix(53, 72) = coeffs(20 - 1);
  coefficient_matrix(54, 14) = coeffs(1 - 1);
  coefficient_matrix(54, 21) = coeffs(12 - 1);
  coefficient_matrix(54, 36) = coeffs(13 - 1);
  coefficient_matrix(54, 38) = coeffs(14 - 1);
  coefficient_matrix(54, 40) = coeffs(15 - 1);
  coefficient_matrix(54, 43) = coeffs(16 - 1);
  coefficient_matrix(54, 54) = coeffs(17 - 1);
  coefficient_matrix(54, 56) = coeffs(18 - 1);
  coefficient_matrix(54, 59) = coeffs(19 - 1);
  coefficient_matrix(54, 71) = coeffs(20 - 1);
  coefficient_matrix(55, 10) = coeffs(1 - 1);
  coefficient_matrix(55, 19) = coeffs(12 - 1);
  coefficient_matrix(55, 32) = coeffs(13 - 1);
  coefficient_matrix(55, 34) = coeffs(14 - 1);
  coefficient_matrix(55, 37) = coeffs(15 - 1);
  coefficient_matrix(55, 51) = coeffs(17 - 1);
  coefficient_matrix(55, 53) = coeffs(18 - 1);
  coefficient_matrix(55, 57) = coeffs(19 - 1);
  coefficient_matrix(55, 69) = coeffs(20 - 1);
  coefficient_matrix(56, 9) = coeffs(1 - 1);
  coefficient_matrix(56, 17) = coeffs(12 - 1);
  coefficient_matrix(56, 30) = coeffs(13 - 1);
  coefficient_matrix(56, 32) = coeffs(14 - 1);
  coefficient_matrix(56, 35) = coeffs(15 - 1);
  coefficient_matrix(56, 41) = coeffs(16 - 1);
  coefficient_matrix(56, 49) = coeffs(17 - 1);
  coefficient_matrix(56, 51) = coeffs(18 - 1);
  coefficient_matrix(56, 55) = coeffs(19 - 1);
  coefficient_matrix(56, 67) = coeffs(20 - 1);
  coefficient_matrix(57, 8) = coeffs(1 - 1);
  coefficient_matrix(57, 16) = coeffs(12 - 1);
  coefficient_matrix(57, 29) = coeffs(13 - 1);
  coefficient_matrix(57, 31) = coeffs(14 - 1);
  coefficient_matrix(57, 40) = coeffs(16 - 1);
  coefficient_matrix(57, 50) = coeffs(18 - 1);
  coefficient_matrix(57, 54) = coeffs(19 - 1);
  coefficient_matrix(58, 6) = coeffs(1 - 1);
  coefficient_matrix(58, 11) = coeffs(12 - 1);
  coefficient_matrix(58, 24) = coeffs(13 - 1);
  coefficient_matrix(58, 26) = coeffs(14 - 1);
  coefficient_matrix(58, 28) = coeffs(15 - 1);
  coefficient_matrix(58, 35) = coeffs(16 - 1);
  coefficient_matrix(58, 44) = coeffs(18 - 1);
  coefficient_matrix(58, 49) = coeffs(19 - 1);
  coefficient_matrix(58, 62) = coeffs(20 - 1);
  coefficient_matrix(59, 17) = coeffs(1 - 1);
  coefficient_matrix(59, 21) = coeffs(21 - 1);
  coefficient_matrix(59, 35) = coeffs(3 - 1);
  coefficient_matrix(59, 40) = coeffs(22 - 1);
  coefficient_matrix(59, 41) = coeffs(23 - 1);
  coefficient_matrix(59, 43) = coeffs(24 - 1);
  coefficient_matrix(59, 54) = coeffs(25 - 1);
  coefficient_matrix(59, 55) = coeffs(26 - 1);
  coefficient_matrix(59, 59) = coeffs(27 - 1);
  coefficient_matrix(59, 71) = coeffs(28 - 1);
  coefficient_matrix(60, 13) = coeffs(1 - 1);
  coefficient_matrix(60, 18) = coeffs(21 - 1);
  coefficient_matrix(60, 30) = coeffs(3 - 1);
  coefficient_matrix(60, 36) = coeffs(22 - 1);
  coefficient_matrix(60, 37) = coeffs(23 - 1);
  coefficient_matrix(60, 42) = coeffs(24 - 1);
  coefficient_matrix(60, 50) = coeffs(25 - 1);
  coefficient_matrix(60, 51) = coeffs(26 - 1);
  coefficient_matrix(60, 56) = coeffs(27 - 1);
  coefficient_matrix(60, 68) = coeffs(28 - 1);
  coefficient_matrix(61, 11) = coeffs(1 - 1);
  coefficient_matrix(61, 16) = coeffs(21 - 1);
  coefficient_matrix(61, 28) = coeffs(3 - 1);
  coefficient_matrix(61, 35) = coeffs(23 - 1);
  coefficient_matrix(61, 40) = coeffs(24 - 1);
  coefficient_matrix(61, 49) = coeffs(26 - 1);
  coefficient_matrix(61, 54) = coeffs(27 - 1);
  coefficient_matrix(62, 7) = coeffs(1 - 1);
  coefficient_matrix(62, 12) = coeffs(21 - 1);
  coefficient_matrix(62, 24) = coeffs(3 - 1);
  coefficient_matrix(62, 29) = coeffs(22 - 1);
  coefficient_matrix(62, 30) = coeffs(23 - 1);
  coefficient_matrix(62, 36) = coeffs(24 - 1);
  coefficient_matrix(62, 44) = coeffs(26 - 1);
  coefficient_matrix(62, 50) = coeffs(27 - 1);
  coefficient_matrix(62, 63) = coeffs(28 - 1);
  coefficient_matrix(63, 20) = coeffs(1 - 1);
  coefficient_matrix(63, 23) = coeffs(29 - 1);
  coefficient_matrix(63, 38) = coeffs(3 - 1);
  coefficient_matrix(63, 42) = coeffs(30 - 1);
  coefficient_matrix(63, 56) = coeffs(33 - 1);
  coefficient_matrix(63, 58) = coeffs(34 - 1);
  coefficient_matrix(63, 61) = coeffs(35 - 1);
  coefficient_matrix(63, 73) = coeffs(36 - 1);
  coefficient_matrix(64, 19) = coeffs(1 - 1);
  coefficient_matrix(64, 22) = coeffs(29 - 1);
  coefficient_matrix(64, 37) = coeffs(3 - 1);
  coefficient_matrix(64, 41) = coeffs(30 - 1);
  coefficient_matrix(64, 55) = coeffs(33 - 1);
  coefficient_matrix(64, 57) = coeffs(34 - 1);
  coefficient_matrix(64, 60) = coeffs(35 - 1);
  coefficient_matrix(64, 72) = coeffs(36 - 1);
  coefficient_matrix(65, 18) = coeffs(1 - 1);
  coefficient_matrix(65, 21) = coeffs(29 - 1);
  coefficient_matrix(65, 36) = coeffs(3 - 1);
  coefficient_matrix(65, 40) = coeffs(30 - 1);
  coefficient_matrix(65, 42) = coeffs(31 - 1);
  coefficient_matrix(65, 43) = coeffs(32 - 1);
  coefficient_matrix(65, 54) = coeffs(33 - 1);
  coefficient_matrix(65, 56) = coeffs(34 - 1);
  coefficient_matrix(65, 59) = coeffs(35 - 1);
  coefficient_matrix(65, 71) = coeffs(36 - 1);
  coefficient_matrix(66, 14) = coeffs(1 - 1);
  coefficient_matrix(66, 18) = coeffs(29 - 1);
  coefficient_matrix(66, 31) = coeffs(3 - 1);
  coefficient_matrix(66, 36) = coeffs(30 - 1);
  coefficient_matrix(66, 38) = coeffs(31 - 1);
  coefficient_matrix(66, 42) = coeffs(32 - 1);
  coefficient_matrix(66, 50) = coeffs(33 - 1);
  coefficient_matrix(66, 52) = coeffs(34 - 1);
  coefficient_matrix(66, 56) = coeffs(35 - 1);
  coefficient_matrix(66, 68) = coeffs(36 - 1);
  coefficient_matrix(67, 12) = coeffs(1 - 1);
  coefficient_matrix(67, 16) = coeffs(29 - 1);
  coefficient_matrix(67, 29) = coeffs(3 - 1);
  coefficient_matrix(67, 36) = coeffs(31 - 1);
  coefficient_matrix(67, 40) = coeffs(32 - 1);
  coefficient_matrix(67, 50) = coeffs(34 - 1);
  coefficient_matrix(67, 54) = coeffs(35 - 1);
  coefficient_matrix(68, 10) = coeffs(1 - 1);
  coefficient_matrix(68, 15) = coeffs(29 - 1);
  coefficient_matrix(68, 27) = coeffs(3 - 1);
  coefficient_matrix(68, 32) = coeffs(30 - 1);
  coefficient_matrix(68, 34) = coeffs(31 - 1);
  coefficient_matrix(68, 39) = coeffs(32 - 1);
  coefficient_matrix(68, 46) = coeffs(33 - 1);
  coefficient_matrix(68, 48) = coeffs(34 - 1);
  coefficient_matrix(68, 53) = coeffs(35 - 1);
  coefficient_matrix(68, 66) = coeffs(36 - 1);
  coefficient_matrix(69, 8) = coeffs(1 - 1);
  coefficient_matrix(69, 12) = coeffs(29 - 1);
  coefficient_matrix(69, 25) = coeffs(3 - 1);
  coefficient_matrix(69, 29) = coeffs(30 - 1);
  coefficient_matrix(69, 31) = coeffs(31 - 1);
  coefficient_matrix(69, 36) = coeffs(32 - 1);
  coefficient_matrix(69, 45) = coeffs(34 - 1);
  coefficient_matrix(69, 50) = coeffs(35 - 1);
  coefficient_matrix(69, 63) = coeffs(36 - 1);
  coefficient_matrix(70, 1) = coeffs(1 - 1);
  coefficient_matrix(70, 5) = coeffs(2 - 1);
  coefficient_matrix(70, 11) = coeffs(3 - 1);
  coefficient_matrix(70, 12) = coeffs(3 - 1);
  coefficient_matrix(70, 13) = coeffs(4 - 1);
  coefficient_matrix(70, 16) = coeffs(5 - 1);
  coefficient_matrix(70, 21) = coeffs(6 - 1);
  coefficient_matrix(70, 35) = coeffs(8 - 1);
  coefficient_matrix(70, 36) = coeffs(9 - 1);
  coefficient_matrix(70, 40) = coeffs(10 - 1);
  coefficient_matrix(70, 54) = coeffs(11 - 1);
  coefficient_matrix(71, 0) = coeffs(1 - 1);
  coefficient_matrix(71, 4) = coeffs(2 - 1);
  coefficient_matrix(71, 7) = coeffs(3 - 1);
  coefficient_matrix(71, 8) = coeffs(3 - 1);
  coefficient_matrix(71, 9) = coeffs(4 - 1);
  coefficient_matrix(71, 12) = coeffs(5 - 1);
  coefficient_matrix(71, 18) = coeffs(6 - 1);
  coefficient_matrix(71, 29) = coeffs(7 - 1);
  coefficient_matrix(71, 30) = coeffs(8 - 1);
  coefficient_matrix(71, 31) = coeffs(9 - 1);
  coefficient_matrix(71, 36) = coeffs(10 - 1);
  coefficient_matrix(71, 50) = coeffs(11 - 1);
  coefficient_matrix(72, 2) = coeffs(1 - 1);
  coefficient_matrix(72, 5) = coeffs(12 - 1);
  coefficient_matrix(72, 12) = coeffs(13 - 1);
  coefficient_matrix(72, 14) = coeffs(14 - 1);
  coefficient_matrix(72, 16) = coeffs(15 - 1);
  coefficient_matrix(72, 21) = coeffs(16 - 1);
  coefficient_matrix(72, 36) = coeffs(18 - 1);
  coefficient_matrix(72, 40) = coeffs(19 - 1);
  coefficient_matrix(72, 54) = coeffs(20 - 1);
  coefficient_matrix(73, 0) = coeffs(1 - 1);
  coefficient_matrix(73, 3) = coeffs(12 - 1);
  coefficient_matrix(73, 7) = coeffs(13 - 1);
  coefficient_matrix(73, 9) = coeffs(14 - 1);
  coefficient_matrix(73, 11) = coeffs(15 - 1);
  coefficient_matrix(73, 17) = coeffs(16 - 1);
  coefficient_matrix(73, 28) = coeffs(17 - 1);
  coefficient_matrix(73, 30) = coeffs(18 - 1);
  coefficient_matrix(73, 35) = coeffs(19 - 1);
  coefficient_matrix(73, 49) = coeffs(20 - 1);
  coefficient_matrix(74, 3) = coeffs(1 - 1);
  coefficient_matrix(74, 5) = coeffs(21 - 1);
  coefficient_matrix(74, 11) = coeffs(3 - 1);
  coefficient_matrix(74, 16) = coeffs(22 - 1);
  coefficient_matrix(74, 17) = coeffs(23 - 1);
  coefficient_matrix(74, 21) = coeffs(24 - 1);
  coefficient_matrix(74, 35) = coeffs(26 - 1);
  coefficient_matrix(74, 40) = coeffs(27 - 1);
  coefficient_matrix(74, 54) = coeffs(28 - 1);
  coefficient_matrix(75, 1) = coeffs(1 - 1);
  coefficient_matrix(75, 4) = coeffs(21 - 1);
  coefficient_matrix(75, 7) = coeffs(3 - 1);
  coefficient_matrix(75, 12) = coeffs(22 - 1);
  coefficient_matrix(75, 13) = coeffs(23 - 1);
  coefficient_matrix(75, 18) = coeffs(24 - 1);
  coefficient_matrix(75, 29) = coeffs(25 - 1);
  coefficient_matrix(75, 30) = coeffs(26 - 1);
  coefficient_matrix(75, 36) = coeffs(27 - 1);
  coefficient_matrix(75, 50) = coeffs(28 - 1);
  coefficient_matrix(76, 4) = coeffs(1 - 1);
  coefficient_matrix(76, 5) = coeffs(29 - 1);
  coefficient_matrix(76, 12) = coeffs(3 - 1);
  coefficient_matrix(76, 16) = coeffs(30 - 1);
  coefficient_matrix(76, 18) = coeffs(31 - 1);
  coefficient_matrix(76, 21) = coeffs(32 - 1);
  coefficient_matrix(76, 36) = coeffs(34 - 1);
  coefficient_matrix(76, 40) = coeffs(35 - 1);
  coefficient_matrix(76, 54) = coeffs(36 - 1);
  coefficient_matrix(77, 2) = coeffs(1 - 1);
  coefficient_matrix(77, 4) = coeffs(29 - 1);
  coefficient_matrix(77, 8) = coeffs(3 - 1);
  coefficient_matrix(77, 12) = coeffs(30 - 1);
  coefficient_matrix(77, 14) = coeffs(31 - 1);
  coefficient_matrix(77, 18) = coeffs(32 - 1);
  coefficient_matrix(77, 29) = coeffs(33 - 1);
  coefficient_matrix(77, 31) = coeffs(34 - 1);
  coefficient_matrix(77, 36) = coeffs(35 - 1);
  coefficient_matrix(77, 50) = coeffs(36 - 1);

  b(0, 1) = coeffs(5 - 1);
  b(0, 4) = coeffs(6 - 1);
  b(0, 5) = coeffs(7 - 1);
  b(0, 6) = coeffs(8 - 1);
  b(0, 7) = coeffs(9 - 1);
  b(0, 8) = coeffs(10 - 1);
  b(0, 9) = coeffs(11 - 1);
  b(1, 0) = coeffs(14 - 1);
  b(1, 1) = coeffs(15 - 1);
  b(1, 4) = coeffs(16 - 1);
  b(1, 5) = coeffs(17 - 1);
  b(1, 7) = coeffs(18 - 1);
  b(1, 8) = coeffs(19 - 1);
  b(1, 9) = coeffs(20 - 1);
  b(2, 1) = coeffs(22 - 1);
  b(2, 2) = coeffs(23 - 1);
  b(2, 4) = coeffs(24 - 1);
  b(2, 5) = coeffs(25 - 1);
  b(2, 6) = coeffs(26 - 1);
  b(2, 8) = coeffs(27 - 1);
  b(2, 9) = coeffs(28 - 1);
  b(3, 1) = coeffs(30 - 1);
  b(3, 3) = coeffs(31 - 1);
  b(3, 4) = coeffs(32 - 1);
  b(3, 5) = coeffs(33 - 1);
  b(3, 7) = coeffs(34 - 1);
  b(3, 8) = coeffs(35 - 1);
  b(3, 9) = coeffs(36 - 1);
  b(4, 1) = coeffs(7 - 1);
  b(4, 2) = coeffs(8 - 1);
  b(4, 3) = coeffs(9 - 1);
  b(4, 4) = coeffs(10 - 1);
  b(4, 8) = coeffs(11 - 1);
  b(5, 0) = coeffs(9 - 1);
  b(5, 3) = coeffs(10 - 1);
  b(5, 7) = coeffs(11 - 1);
  b(6, 1) = coeffs(10 - 1);
  b(6, 5) = coeffs(11 - 1);
  b(7, 1) = coeffs(17 - 1);
  b(7, 3) = coeffs(18 - 1);
  b(7, 4) = coeffs(19 - 1);
  b(7, 8) = coeffs(20 - 1);
  b(8, 2) = coeffs(19 - 1);
  b(8, 6) = coeffs(20 - 1);
  b(9, 1) = coeffs(19 - 1);
  b(9, 5) = coeffs(20 - 1);
  b(10, 1) = coeffs(25 - 1);
  b(10, 2) = coeffs(26 - 1);
  b(10, 4) = coeffs(27 - 1);
  b(10, 8) = coeffs(28 - 1);
  b(11, 3) = coeffs(27 - 1);
  b(11, 7) = coeffs(28 - 1);
  b(12, 1) = coeffs(27 - 1);
  b(12, 5) = coeffs(28 - 1);
  b(13, 1) = coeffs(33 - 1);
  b(13, 3) = coeffs(34 - 1);
  b(13, 4) = coeffs(35 - 1);
  b(13, 8) = coeffs(36 - 1);
  b(14, 0) = coeffs(34 - 1);
  b(14, 3) = coeffs(35 - 1);
  b(14, 7) = coeffs(36 - 1);
  b(15, 2) = coeffs(35 - 1);
  b(15, 6) = coeffs(36 - 1);
  b(16, 1) = coeffs(35 - 1);
  b(16, 5) = coeffs(36 - 1);
  b(17, 4) = coeffs(11 - 1);
  b(18, 3) = coeffs(11 - 1);
  b(19, 1) = coeffs(11 - 1);
  b(20, 0) = coeffs(11 - 1);
  b(23, 4) = coeffs(20 - 1);
  b(24, 3) = coeffs(20 - 1);
  b(25, 2) = coeffs(20 - 1);
  b(26, 1) = coeffs(20 - 1);
  b(30, 4) = coeffs(28 - 1);
  b(31, 3) = coeffs(28 - 1);
  b(32, 1) = coeffs(28 - 1);
  b(33, 0) = coeffs(28 - 1);
  b(36, 4) = coeffs(36 - 1);
  b(37, 3) = coeffs(36 - 1);
  b(38, 2) = coeffs(36 - 1);
  b(39, 1) = coeffs(36 - 1);
  b(40, 0) = coeffs(36 - 1);

  const Matrix<double, 78, 10> solved_mat =
      coefficient_matrix.partialPivLu().solve(b);

  Matrix<double, 10, 10> eigen_mat = Matrix<double, 10, 10>::Zero();
  eigen_mat(0, 1) = 1.0;
  eigen_mat(1, 5) = 1.0;
  eigen_mat(2, 6) = 1.0;
  eigen_mat(3, 7) = 1.0;
  eigen_mat(4, 8) = 1.0;
  eigen_mat.row(5) = -solved_mat.row(74).reverse();
  eigen_mat.row(6) = -solved_mat.row(73).reverse();
  eigen_mat.row(7) = -solved_mat.row(72).reverse();
  eigen_mat.row(8) = -solved_mat.row(71).reverse();
  eigen_mat.row(9) = -solved_mat.row(70).reverse();

  Eigen::EigenSolver<Matrix<double, 10, 10> > eigensolver(eigen_mat);
  const Eigen::EigenSolver<Matrix<double, 10, 10> >::EigenvectorsType
      eigenvectors = eigensolver.eigenvectors();

  Eigen::EigenSolver<Eigen::Matrix<double, 4, 10> >::EigenvectorsType sol;
  sol.row(0).array() =
      eigenvectors.row(1).array() / eigenvectors.row(0).array();
  sol.row(1).array() =
      eigenvectors.row(2).array() / eigenvectors.row(0).array();
  sol.row(2).array() =
      eigenvectors.row(3).array() / eigenvectors.row(0).array();
  sol.row(3).array() =
      eigenvectors.row(4).array() / eigenvectors.row(0).array();

  for (int i = 0; i < 10; i++) {
    if (sol(3, i).imag() != 0 || sol(3, i).real() < 0 ||
        std::isnan(sol(3, i).real())) {
      continue;
    }

    f->push_back(sqrt(sol(3, i).real()));
    depths->push_back(sol.col(i).real().head<3>().reverse());
  }
  return true;
}

}  // namespace theia
