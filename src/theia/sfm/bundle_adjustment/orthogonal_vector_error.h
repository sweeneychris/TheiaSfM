// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_BUNDLE_ADJUSTMENT_ORTHOGONAL_VECTOR_ERROR_H_
#define THEIA_SFM_BUNDLE_ADJUSTMENT_ORTHOGONAL_VECTOR_ERROR_H_

#include <ceres/ceres.h>

#include <Eigen/Core>

namespace theia {

struct OrthogonalVectorError {
 public:
  explicit OrthogonalVectorError(const Eigen::Vector3d& orthogonal_vector)
      : orthogonal_vector_(orthogonal_vector) {}

  template <typename T>
  bool operator()(const T* vector, T* residual) const {
    residual[0] = T(orthogonal_vector_[0]) * vector[0] +
                  T(orthogonal_vector_[1]) * vector[1] +
                  T(orthogonal_vector_[2]) * vector[2];
    return true;
  }

  static ceres::CostFunction* Create(const Eigen::Vector3d& orthogonal_vector) {
    static const int kParameterSize = 3;
    return new ceres::AutoDiffCostFunction<OrthogonalVectorError,
                                           1,
                                           kParameterSize>(
        new OrthogonalVectorError(orthogonal_vector));
  }

 private:
  const Eigen::Vector3d orthogonal_vector_;
};

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_ORTHOGONAL_VECTOR_ERROR_H_
