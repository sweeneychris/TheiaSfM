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

#include "theia/sfm/pose/test_util.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <glog/logging.h>
#include "theia/util/random.h"

namespace theia {

using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// Adds noise to the 3D point passed in.
void AddNoiseToPoint(const double noise_factor, Vector3d* point) {
  *point += Vector3d((-0.5 + RandDouble(0.0, 1.0)) * 2.0 * noise_factor,
                     (-0.5 + RandDouble(0.0, 1.0)) * 2.0 * noise_factor,
                     (-0.5 + RandDouble(0.0, 1.0)) * 2.0 * noise_factor);
}

// Adds noise to the ray i.e. the projection of the point.
void AddNoiseToProjection(const double noise_factor, Vector2d* ray) {
  *ray += noise_factor * Vector2d::Random();
}

// Basically, transform the point to be at (0, 0, -1), add noise, then reverse
// the transformation.
void AddNoiseToRay(const double std_dev, Vector3d* proj) {
  const double scale = proj->norm();
  const double noise_x = (-0.5 + theia::RandDouble(0.0, 1.0)) * 2.0 * std_dev;
  const double noise_y = (-0.5 + theia::RandDouble(0.0, 1.0)) * 2.0 * std_dev;

  Eigen::Quaterniond rot =
      Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(0, 0, 1.0), *proj);
  Eigen::Vector3d noisy_point(noise_x, noise_y, 1);
  noisy_point *= scale;

  *proj = (rot * noisy_point).normalized();
}

void AddGaussianNoise(const double noise_factor, Vector2d* ray) {
  const double noise_x = RandGaussian(0.0, noise_factor);
  const double noise_y = RandGaussian(0.0, noise_factor);
  *ray = Vector2d(ray->x() + noise_x, ray->y() + noise_y);
}

void CreateRandomPointsInFrustum(const double near_plane_width,
                                 const double near_plane_height,
                                 const double near_plane_depth,
                                 const double far_plane_depth,
                                 const int num_points,
                                 std::vector<Eigen::Vector3d>* random_points) {
  random_points->reserve(num_points);
  for (int i = 0; i < num_points; i++) {
    const double rand_depth = RandDouble(near_plane_depth, far_plane_depth);
    const double x_radius = near_plane_width * rand_depth / near_plane_depth;
    const double y_radius = near_plane_height * rand_depth / near_plane_depth;
    Vector3d rand_point(RandDouble(-x_radius, x_radius),
                        RandDouble(-y_radius, y_radius), rand_depth);
    random_points->push_back(rand_point);
  }
}

}  // namespace theia
