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

#ifndef THEIA_SFM_BUNDLE_ADJUSTMENT_UNIT_NORM_THREE_VECTOR_PARAMETERIZATION_H_
#define THEIA_SFM_BUNDLE_ADJUSTMENT_UNIT_NORM_THREE_VECTOR_PARAMETERIZATION_H_

#include <ceres/ceres.h>

namespace theia {

// A parameterization for Ceres that keeps a 3-dimensional vector as a unit-norm
// vector throughout optimization.
struct UnitNormThreeVectorParameterization {
  template<typename T>
  bool operator()(const T* x, const T* delta, T* x_plus_delta) const {
    x_plus_delta[0] = x[0] + delta[0];
    x_plus_delta[1] = x[1] + delta[1];
    x_plus_delta[2] = x[2] + delta[2];
    const T sq_norm = x_plus_delta[0] * x_plus_delta[0] +
                      x_plus_delta[1] * x_plus_delta[1] +
                      x_plus_delta[2] * x_plus_delta[2];
    if (sq_norm > T(0.0)) {
      const T norm = sqrt(sq_norm);
      x_plus_delta[0] /= norm;
      x_plus_delta[1] /= norm;
      x_plus_delta[2] /= norm;
    }

    return true;
  }
};

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_UNIT_NORM_THREE_VECTOR_PARAMETERIZATION_H_
