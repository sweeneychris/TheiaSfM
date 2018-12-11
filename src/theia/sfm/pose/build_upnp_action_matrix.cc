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

#include "theia/sfm/pose/build_upnp_action_matrix.h"

#include <Eigen/Core>
#include <vector>

namespace theia {
typedef Eigen::Matrix<double, 8, 8> Matrix8d;
typedef Eigen::Matrix<double, 10, 10> Matrix10d;
typedef Eigen::Matrix<double, 16, 16> Matrix16d;
typedef Eigen::Matrix<double, 10, 1> Vector10d;

// TODO(vfragoso): Document me!
// Implementation based on:
// OpenGV file: src/absolute_pose/modules/upnp4.cpp
Matrix8d BuildActionMatrixUsingSymmetry(const Matrix10d& a_matrix,
                                        const Vector10d& b_vector,
                                        const double gamma) {
  Matrix8d action_matrix;
  const Matrix10d& M = a_matrix;
  const Vector10d& C = b_vector;
  Eigen::Matrix<double, 8, 24> M1 = Eigen::Matrix<double, 8, 24>::Zero();
  M1(0,0) = 4*M(0,0); M1(0,1) = 4*M(4,0)+2*M(0,4); M1(0,2) = 2*M(4,4)+4*M(1,0); M1(0,3) = 2*M(1,4); M1(0,4) = 4*M(5,0)+2*M(0,5); M1(0,5) = 4*M(7,0)+2*M(5,4)+2*M(4,5); M1(0,6) = 2*M(7,4)+2*M(1,5); M1(0,7) = 2*M(5,5)+4*M(2,0); M1(0,8) = 2*M(7,5)+2*M(2,4); M1(0,9) = 2*M(2,5); M1(0,10) = 4*M(6,0)+2*M(0,6); M1(0,11) = 4*M(8,0)+2*M(6,4)+2*M(4,6); M1(0,12) = 2*M(8,4)+2*M(1,6); M1(0,13) = 4*M(9,0)+2*M(6,5)+2*M(5,6); M1(0,14) = 2*M(9,4)+2*M(8,5)+2*M(7,6); M1(0,15) = 2*M(9,5)+2*M(2,6); M1(0,16) = 2*M(6,6)+4*M(3,0); M1(0,17) = 2*M(8,6)+2*M(3,4); M1(0,18) = 2*M(9,6)+2*M(3,5); M1(0,19) = 2*M(3,6); M1(0,20) = 4*C(0,0); M1(0,21) = 2*C(0,4); M1(0,22) = 2*C(0,5); M1(0,23) = 2*C(0,6); 
  M1(1,0) = 2*M(0,4); M1(1,1) = 2*M(4,4)+4*M(0,1); M1(1,2) = 4*M(4,1)+2*M(1,4); M1(1,3) = 4*M(1,1); M1(1,4) = 2*M(5,4)+2*M(0,7); M1(1,5) = 2*M(7,4)+4*M(5,1)+2*M(4,7); M1(1,6) = 4*M(7,1)+2*M(1,7); M1(1,7) = 2*M(5,7)+2*M(2,4); M1(1,8) = 2*M(7,7)+4*M(2,1); M1(1,9) = 2*M(2,7); M1(1,10) = 2*M(6,4)+2*M(0,8); M1(1,11) = 2*M(8,4)+4*M(6,1)+2*M(4,8); M1(1,12) = 4*M(8,1)+2*M(1,8); M1(1,13) = 2*M(9,4)+2*M(6,7)+2*M(5,8); M1(1,14) = 4*M(9,1)+2*M(8,7)+2*M(7,8); M1(1,15) = 2*M(9,7)+2*M(2,8); M1(1,16) = 2*M(6,8)+2*M(3,4); M1(1,17) = 2*M(8,8)+4*M(3,1); M1(1,18) = 2*M(9,8)+2*M(3,7); M1(1,19) = 2*M(3,8); M1(1,20) = 2*C(0,4); M1(1,21) = 4*C(0,1); M1(1,22) = 2*C(0,7); M1(1,23) = 2*C(0,8); 
  M1(2,0) = 2*M(0,5); M1(2,1) = 2*M(4,5)+2*M(0,7); M1(2,2) = 2*M(4,7)+2*M(1,5); M1(2,3) = 2*M(1,7); M1(2,4) = 2*M(5,5)+4*M(0,2); M1(2,5) = 2*M(7,5)+2*M(5,7)+4*M(4,2); M1(2,6) = 2*M(7,7)+4*M(1,2); M1(2,7) = 4*M(5,2)+2*M(2,5); M1(2,8) = 4*M(7,2)+2*M(2,7); M1(2,9) = 4*M(2,2); M1(2,10) = 2*M(6,5)+2*M(0,9); M1(2,11) = 2*M(8,5)+2*M(6,7)+2*M(4,9); M1(2,12) = 2*M(8,7)+2*M(1,9); M1(2,13) = 2*M(9,5)+4*M(6,2)+2*M(5,9); M1(2,14) = 2*M(9,7)+4*M(8,2)+2*M(7,9); M1(2,15) = 4*M(9,2)+2*M(2,9); M1(2,16) = 2*M(6,9)+2*M(3,5); M1(2,17) = 2*M(8,9)+2*M(3,7); M1(2,18) = 2*M(9,9)+4*M(3,2); M1(2,19) = 2*M(3,9); M1(2,20) = 2*C(0,5); M1(2,21) = 2*C(0,7); M1(2,22) = 4*C(0,2); M1(2,23) = 2*C(0,9); 
  M1(3,0) = 2*M(0,6); M1(3,1) = 2*M(4,6)+2*M(0,8); M1(3,2) = 2*M(4,8)+2*M(1,6); M1(3,3) = 2*M(1,8); M1(3,4) = 2*M(5,6)+2*M(0,9); M1(3,5) = 2*M(7,6)+2*M(5,8)+2*M(4,9); M1(3,6) = 2*M(7,8)+2*M(1,9); M1(3,7) = 2*M(5,9)+2*M(2,6); M1(3,8) = 2*M(7,9)+2*M(2,8); M1(3,9) = 2*M(2,9); M1(3,10) = 2*M(6,6)+4*M(0,3); M1(3,11) = 2*M(8,6)+2*M(6,8)+4*M(4,3); M1(3,12) = 2*M(8,8)+4*M(1,3); M1(3,13) = 2*M(9,6)+2*M(6,9)+4*M(5,3); M1(3,14) = 2*M(9,8)+2*M(8,9)+4*M(7,3); M1(3,15) = 2*M(9,9)+4*M(2,3); M1(3,16) = 4*M(6,3)+2*M(3,6); M1(3,17) = 4*M(8,3)+2*M(3,8); M1(3,18) = 4*M(9,3)+2*M(3,9); M1(3,19) = 4*M(3,3); M1(3,20) = 2*C(0,6); M1(3,21) = 2*C(0,8); M1(3,22) = 2*C(0,9); M1(3,23) = 4*C(0,3); 
  M1(4,0) = 1; M1(4,2) = 1; M1(4,7) = 1; M1(4,16) = 1; M1(4,20) = -1; 
  M1(5,1) = 1; M1(5,3) = 1; M1(5,8) = 1; M1(5,17) = 1; M1(5,21) = -1; 
  M1(6,4) = 1; M1(6,6) = 1; M1(6,9) = 1; M1(6,18) = 1; M1(6,22) = -1; 
  M1(7,10) = 1; M1(7,12) = 1; M1(7,15) = 1; M1(7,19) = 1; M1(7,23) = -1; 
  
  const Eigen::Matrix<double, 7, 7> M1_part1 = M1.block<7, 7>(0, 0);
  const Eigen::Matrix<double, 7, 24> M1_part2 =
      M1_part1.inverse() * M1.block<7, 24>(0, 0);
  M1.block<7, 24>(0, 0) = M1_part2;

  // Some more cancellation in column 10.
  for (int i = 6; i >= 0; --i) {
    M1.row(i) -= (M1(i, 10) / M1(7, 10)) * M1.row(7);
  }
  
  return action_matrix;
}

// TODO(vfragoso): Document me!
// Implementation based on:
// OpenGV file: src/absolute_pose/modules/upnp2.cpp
Matrix16d BuildActionMatrix(const Matrix10d& a_matrix,
                            const Vector10d& b_vector,
                            const double gamma) {
  Matrix16d action_matrix;
  const Matrix10d& M = a_matrix;
  const Vector10d& C = b_vector;
  Eigen::Matrix<double, 5, 29> M1 = Eigen::Matrix<double, 5, 29>::Zero();
  M1(0,0) = 4*M(0,0); M1(0,1) = 4*M(4,0)+2*M(0,4); M1(0,2) = 2*M(4,4)+4*M(1,0); M1(0,3) = 2*M(1,4); M1(0,4) = 4*M(5,0)+2*M(0,5); M1(0,5) = 4*M(7,0)+2*M(5,4)+2*M(4,5); M1(0,6) = 2*M(7,4)+2*M(1,5); M1(0,7) = 2*M(5,5)+4*M(2,0); M1(0,8) = 2*M(7,5)+2*M(2,4); M1(0,9) = 2*M(2,5); M1(0,10) = 4*M(6,0)+2*M(0,6); M1(0,11) = 4*M(8,0)+2*M(6,4)+2*M(4,6); M1(0,12) = 2*M(8,4)+2*M(1,6); M1(0,13) = 4*M(9,0)+2*M(6,5)+2*M(5,6); M1(0,14) = 2*M(9,4)+2*M(8,5)+2*M(7,6); M1(0,15) = 2*M(9,5)+2*M(2,6); M1(0,16) = 2*M(6,6)+4*M(3,0); M1(0,17) = 2*M(8,6)+2*M(3,4); M1(0,18) = 2*M(9,6)+2*M(3,5); M1(0,19) = 2*M(3,6); M1(0,24) = 4*C(0,0); M1(0,25) = 2*C(0,4); M1(0,26) = 2*C(0,5); M1(0,27) = 2*C(0,6); 
  M1(1,0) = 2*M(0,4); M1(1,1) = 2*M(4,4)+4*M(0,1); M1(1,2) = 4*M(4,1)+2*M(1,4); M1(1,3) = 4*M(1,1); M1(1,4) = 2*M(5,4)+2*M(0,7); M1(1,5) = 2*M(7,4)+4*M(5,1)+2*M(4,7); M1(1,6) = 4*M(7,1)+2*M(1,7); M1(1,7) = 2*M(5,7)+2*M(2,4); M1(1,8) = 2*M(7,7)+4*M(2,1); M1(1,9) = 2*M(2,7); M1(1,10) = 2*M(6,4)+2*M(0,8); M1(1,11) = 2*M(8,4)+4*M(6,1)+2*M(4,8); M1(1,12) = 4*M(8,1)+2*M(1,8); M1(1,13) = 2*M(9,4)+2*M(6,7)+2*M(5,8); M1(1,14) = 4*M(9,1)+2*M(8,7)+2*M(7,8); M1(1,15) = 2*M(9,7)+2*M(2,8); M1(1,16) = 2*M(6,8)+2*M(3,4); M1(1,17) = 2*M(8,8)+4*M(3,1); M1(1,18) = 2*M(9,8)+2*M(3,7); M1(1,19) = 2*M(3,8); M1(1,24) = 2*C(0,4); M1(1,25) = 4*C(0,1); M1(1,26) = 2*C(0,7); M1(1,27) = 2*C(0,8); 
  M1(2,0) = 2*M(0,5); M1(2,1) = 2*M(4,5)+2*M(0,7); M1(2,2) = 2*M(4,7)+2*M(1,5); M1(2,3) = 2*M(1,7); M1(2,4) = 2*M(5,5)+4*M(0,2); M1(2,5) = 2*M(7,5)+2*M(5,7)+4*M(4,2); M1(2,6) = 2*M(7,7)+4*M(1,2); M1(2,7) = 4*M(5,2)+2*M(2,5); M1(2,8) = 4*M(7,2)+2*M(2,7); M1(2,9) = 4*M(2,2); M1(2,10) = 2*M(6,5)+2*M(0,9); M1(2,11) = 2*M(8,5)+2*M(6,7)+2*M(4,9); M1(2,12) = 2*M(8,7)+2*M(1,9); M1(2,13) = 2*M(9,5)+4*M(6,2)+2*M(5,9); M1(2,14) = 2*M(9,7)+4*M(8,2)+2*M(7,9); M1(2,15) = 4*M(9,2)+2*M(2,9); M1(2,16) = 2*M(6,9)+2*M(3,5); M1(2,17) = 2*M(8,9)+2*M(3,7); M1(2,18) = 2*M(9,9)+4*M(3,2); M1(2,19) = 2*M(3,9); M1(2,24) = 2*C(0,5); M1(2,25) = 2*C(0,7); M1(2,26) = 4*C(0,2); M1(2,27) = 2*C(0,9); 
  M1(3,0) = 2*M(0,6); M1(3,1) = 2*M(4,6)+2*M(0,8); M1(3,2) = 2*M(4,8)+2*M(1,6); M1(3,3) = 2*M(1,8); M1(3,4) = 2*M(5,6)+2*M(0,9); M1(3,5) = 2*M(7,6)+2*M(5,8)+2*M(4,9); M1(3,6) = 2*M(7,8)+2*M(1,9); M1(3,7) = 2*M(5,9)+2*M(2,6); M1(3,8) = 2*M(7,9)+2*M(2,8); M1(3,9) = 2*M(2,9); M1(3,10) = 2*M(6,6)+4*M(0,3); M1(3,11) = 2*M(8,6)+2*M(6,8)+4*M(4,3); M1(3,12) = 2*M(8,8)+4*M(1,3); M1(3,13) = 2*M(9,6)+2*M(6,9)+4*M(5,3); M1(3,14) = 2*M(9,8)+2*M(8,9)+4*M(7,3); M1(3,15) = 2*M(9,9)+4*M(2,3); M1(3,16) = 4*M(6,3)+2*M(3,6); M1(3,17) = 4*M(8,3)+2*M(3,8); M1(3,18) = 4*M(9,3)+2*M(3,9); M1(3,19) = 4*M(3,3); M1(3,24) = 2*C(0,6); M1(3,25) = 2*C(0,8); M1(3,26) = 2*C(0,9); M1(3,27) = 4*C(0,3); 
  M1(4,20) = 1; M1(4,21) = 1; M1(4,22) = 1; M1(4,23) = 1; M1(4,28) = -1;

  const Eigen::Matrix4d M1_part1 = M1.block<4, 4>(0, 0);
  const Eigen::Matrix<double, 4, 29> M1_part2 =
      M1_part1.inverse() * M1.block<4, 29>(0, 0);
  M1.block<4, 29>(0, 0) = M1_part2;

  // TODO(vfragoso): Investigate what GaussJordan code is doing and replace it
  // with Eigen code if possible.
  return action_matrix;
}

}  // namespace theia
