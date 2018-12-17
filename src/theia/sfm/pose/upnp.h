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

#ifndef THEIA_SFM_POSE_UPNP_H_
#define THEIA_SFM_POSE_UPNP_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>

namespace theia {

// TODO(vfragoso): Document me!
struct UpnpCostParameters {
  UpnpCostParameters() {
    a_matrix.setZero();
    b_vector.setZero();
    gamma = 0.0;
  }
  ~UpnpCostParameters() = default;

  Eigen::Matrix<double, 10, 10> a_matrix;
  Eigen::Matrix<double, 10, 1> b_vector;
  double gamma;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// TODO(vfragoso): Document me!
double EvaluateUpnpCost(const UpnpCostParameters& parameters,
                        const Eigen::Quaterniond& rotation);

// TODO(vfragoso): Document me!
UpnpCostParameters ComputeUpnpCostParameters(
    const std::vector<Eigen::Vector3d>& ray_origins,
    const std::vector<Eigen::Vector3d>& ray_directions,
    const std::vector<Eigen::Vector3d>& world_points);

// TODO(vfragoso): Document me!
UpnpCostParameters Upnp(const std::vector<Eigen::Vector3d>& ray_origins,
                        const std::vector<Eigen::Vector3d>& ray_directions,
                        const std::vector<Eigen::Vector3d>& world_points,
                        std::vector<Eigen::Quaterniond>* solution_rotations,
                        std::vector<Eigen::Vector3d>* solution_translations);

}  // namespace theia

#endif  // THEIA_SFM_POSE_UPNP_H_
