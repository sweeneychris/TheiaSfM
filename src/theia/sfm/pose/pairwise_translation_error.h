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

#ifndef THEIA_SFM_POSE_PAIRWISE_TRANSLATION_ERROR_H_
#define THEIA_SFM_POSE_PAIRWISE_TRANSLATION_ERROR_H_

#include <ceres/ceres.h>
#include <Eigen/Core>

namespace theia {

// Computes the error between a translation direction and the direction formed
// from two positions such that (c_j - c_i) - scalar * t_ij is minimized.
struct PairwiseTranslationError {
  PairwiseTranslationError(const Eigen::Vector3d& translation_direction,
                           const double weight);

  // The error is given by the position error described above.
  template <typename T>
  bool operator()(const T* position1, const T* position2, const T* scalar,
                  T* residuals) const;

  static ceres::CostFunction* Create(
      const Eigen::Vector3d& translation_direction, const double weight);

  const Eigen::Vector3d translation_direction_;
  const double weight_;
};

template <typename T>
bool PairwiseTranslationError::operator() (const T* position1,
                                           const T* position2,
                                           const T* scalar,
                                           T* residuals) const {
  residuals[0] = T(weight_) * (position2[0] - position1[0] -
                               scalar[0] * T(translation_direction_[0]));
  residuals[1] = T(weight_) * (position2[1] - position1[1] -
                               scalar[0] * T(translation_direction_[1]));
  residuals[2] = T(weight_) * (position2[2] - position1[2] -
                               scalar[0] * T(translation_direction_[2]));
  return true;
}

}  // namespace theia

#endif  // THEIA_SFM_POSE_PAIRWISE_TRANSLATION_ERROR_H_
