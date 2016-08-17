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

#ifndef THEIA_SFM_POSE_TEST_UTIL_H_
#define THEIA_SFM_POSE_TEST_UTIL_H_

#include <Eigen/Core>
#include <vector>

#include "theia/util/random.h"

namespace theia {

// Adds noise to the 3D point passed in.
void AddNoiseToPoint(const double noise_factor,
                     RandomNumberGenerator* rng,
                     Eigen::Vector3d* point);

// Adds noise to the ray i.e. the projection of the point.
void AddNoiseToProjection(const double noise_factor,
                          RandomNumberGenerator* rng,
                          Eigen::Vector2d* point);

// Adds noise to the image ray.
void AddNoiseToRay(const double std_dev,
                   RandomNumberGenerator* rng,
                   Eigen::Vector3d* proj);

void AddGaussianNoise(const double noise_factor,
                      RandomNumberGenerator* rng,
                      Eigen::Vector3d* ray);

// Creates points that are randomly distributed within a viewing frustum.
void CreateRandomPointsInFrustum(const double near_plane_width,
                                 const double near_plane_height,
                                 const double near_plane_depth,
                                 const double far_plane_depth,
                                 const int num_points,
                                 RandomNumberGenerator* rng,
                                 std::vector<Eigen::Vector3d>* random_points);
}  // namespace theia

#endif  // THEIA_SFM_POSE_TEST_UTIL_H_
