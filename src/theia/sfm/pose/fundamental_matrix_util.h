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

#ifndef THEIA_SFM_POSE_FUNDAMENTAL_MATRIX_UTIL_H_
#define THEIA_SFM_POSE_FUNDAMENTAL_MATRIX_UTIL_H_

namespace theia {

// Given a fundmental matrix, decompose the fundmental matrix and recover focal
// lengths f1, f2 >0. This assumes a principal point of (0, 0) for both cameras.
//
// Returns true on success, false otherwise.
bool FocalLengthsFromFundamentalMatrix(const double fmatrix[3 * 3],
                                       double* focal_length1,
                                       double* focal_length2);

// Given a fundamental matrix that relates two cameras with the same intrinsics,
// extract the shared focal length. This assumes that the fundamental matrix was
// computed with the effect of all non-focal length intrinsics (e.g., principal
// point, aspect ratio, etc.) removed.
bool SharedFocalLengthsFromFundamentalMatrix(const double fmatrix[3 * 3],
                                             double* focal_length);

// Computes the projection matrices corresponding to the fundamental matrix such
// that if y^t * F * x = 0, pmatrix1 corresponds to the camera observing y and
// pmatrix2 corresponds to the camera observing x.
void ProjectionMatricesFromFundamentalMatrix(const double fmatrix[3 * 3],
                                             double pmatrix1[3 * 4],
                                             double pmatrix2[3 * 4]);

// Constructs projection matrices from the input fundamental matrix. The
// fundamental matrix is such that point1^t * fmatrix * point2 = 0 for point1 in
// the image corresponding to pmatrix1 and point2 in the image corresponding to
// pmatrix2.
void FundamentalMatrixFromProjectionMatrices(const double pmatrix1[3 * 4],
                                             const double pmatrix2[3 * 4],
                                             double fmatrix[3 * 3]);

// Extracts the essential matrix such that diag([f2 f2 1]) F diag[f1 f1 1]) is a
// valid essential matrix.
void EssentialMatrixFromFundamentalMatrix(const double fmatrix[3 * 3],
                                          const double focal_length1,
                                          const double focal_length2,
                                          double ematrix[3 * 3]);

// Composes a fundamental matrix such that:
//    F = K_2^-1 * [t]_x * R * K_2^-1
// where K_i is the calibration matrix for image i. The fundamental matrix F is
// thus the matrix that transfers points from image 1 to lines in image 2.
void ComposeFundamentalMatrix(const double focal_length1,
                              const double focal_length2,
                              const double rotation[3 * 3],
                              const double translation[3],
                              double fmatrix[3 * 3]);

}  // namespace theia

#endif  // THEIA_SFM_POSE_FUNDAMENTAL_MATRIX_UTIL_H_
