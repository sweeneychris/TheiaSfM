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

#ifndef THEIA_SFM_TRANSFORMATION_ALIGN_POINT_CLOUDS_H_
#define THEIA_SFM_TRANSFORMATION_ALIGN_POINT_CLOUDS_H_

#include <Eigen/Core>
#include <vector>

namespace theia {

// Computes the orientation, position, and scale factor for the transformation
// between two corresponding 3D point sets A and B such as they are related by:
//
//     B = s * R * A + t
//
// where A is "left" and B is "right". Implementation is based on the paper by
// Umeyama "Least-squares estimation of transformation parameters between two
// point patterns".
void AlignPointCloudsUmeyama(const std::vector<Eigen::Vector3d>& left,
                             const std::vector<Eigen::Vector3d>& right,
                             Eigen::Matrix3d* rotation,
                             Eigen::Vector3d* translation, double* scale);

// Adds weights for each match
// The previous objective function E = Sum(|| yi - (c.R.xi + T) ||^2) becomes
// Sum(wi * || yi - (c.R.xi + T) ||^2)
// The weights must be positive
void AlignPointCloudsUmeyamaWithWeights(
    const std::vector<Eigen::Vector3d>& left,
    const std::vector<Eigen::Vector3d>& right,
    const std::vector<double>& weights, Eigen::Matrix3d* rotation,
    Eigen::Vector3d* translation, double* scale);

}  // namespace theia

#endif  // THEIA_SFM_TRANSFORMATION_ALIGN_POINT_CLOUDS_H_
