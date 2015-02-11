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

#ifndef THEIA_SFM_POSE_DLS_IMPL_H_
#define THEIA_SFM_POSE_DLS_IMPL_H_

#include <Eigen/Core>

namespace theia {
using Eigen::Matrix;
using Eigen::Vector3d;

// Put these methods in a nested namespace so that they are not part of the
// common public theia namespace.
namespace dls_impl {

// Transforms a 3-vector in a 3x9 matrix such that:
// R * v = LeftMultiplyMatrix(v) * vec(R)
// Where R is a rotation matrix and vec(R) converts R to a 9x1 vector.
Eigen::Matrix<double, 3, 9> LeftMultiplyMatrix(const Eigen::Vector3d& v);

// Extracts the coefficients of the Jacobians of the LS cost function (which is
// parameterized by the 3 rotation coefficients s1, s2, s3).
void ExtractJacobianCoefficients(
    const Eigen::Matrix<double, 9, 9>& ls_cost_coefficients,
    double f1_coeff[20], double f2_coeff[20], double f3_coeff[20]);

// Constructs a Macaulay matrix to solve the system of equations using the
// polynomial coefficients from the jacobians.
Eigen::Matrix<double, 120, 120> CreateMacaulayMatrix(
    const double f1_coeff[20], const double f2_coeff[20],
    const double f3_coeff[20], const double rand_term[4]);

}  // namespace dls_impl
}  // namespace theia

#endif  // THEIA_SFM_POSE_DLS_IMPL_H_
