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
#include "theia/alignment/alignment.h"

// TODO(vfragoso): Document me!
namespace theia {

class Upnp {
  // Useful aliases for data types.
  using Matrix10d = Eigen::Matrix<double, 10, 10>;
  using Vector10d = Eigen::Matrix<double, 10, 1>;
  using RowMajorMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  // Dimensions for the minimal-sample template matrix.
  const int kMinimalSampleTemplateRows = 395;
  const int kMinimalSampleTemplateCols = 412;

  // Dimensions for the non-minimal-sample template matrix.
  const int kNonMinimalSampleTemplateRows = 141;
  const int kNonMinimalSampleTemplateCols = 149;

 public:
  // Constructors and destructors.
  Upnp() = default;
  ~Upnp() = default;

  // TODO(vfragoso): Document me!
  struct CostParameters {
    CostParameters() {
      quadratic_penalty_mat.setZero();
      linear_penalty_vector.setZero();
      gamma = 0.0;
    }
    ~CostParameters() = default;

    Matrix10d quadratic_penalty_mat;
    Vector10d linear_penalty_vector;
    double gamma;
  };

  // Estimates poses.
  // TODO(vfragoso): Document me!
  bool EstimatePose(const std::vector<Eigen::Vector3d>& ray_origins,
                    const std::vector<Eigen::Vector3d>& ray_directions,
                    const std::vector<Eigen::Vector3d>& world_points,
                    std::vector<Eigen::Quaterniond>* solution_rotations,
                    std::vector<Eigen::Vector3d>* solution_translations);

  // Getters.
  const CostParameters& cost_params() const {
    return cost_params_;
  }

  // Helper functions.
  // Evaluates the Upnp cost function given a rotation and the cost parameters.
  static double EvaluateCost(const CostParameters& parameters,
                             const Eigen::Quaterniond& rotation);

 private:
  // Cost parameters.
  CostParameters cost_params_;

  // Template matrices.
  RowMajorMatrixXd minimal_sample_template_matrix_;
  RowMajorMatrixXd non_minimal_sample_template_matrix_;

  // Computes the entries of the cost parameters given 2D-3D correspondences.
  std::vector<Eigen::Matrix3d> ComputeCostParameters(
    const std::vector<Eigen::Vector3d>& ray_origins,
    const std::vector<Eigen::Vector3d>& ray_directions,
    const std::vector<Eigen::Vector3d>& world_points);

  // Compute rotations.
  std::vector<Eigen::Quaterniond> ComputeRotations(
      const int num_correspondences);

  // Solve Upnp polynomial solver from minimal sample.
  std::vector<Eigen::Quaterniond> SolveForRotationsFromMinimalSample();

  // Solve Upnp polynomial solver from a non-minimal sample.
  std::vector<Eigen::Quaterniond> SolveForRotationsFromNonMinimalSample();
};

// Estimates the pose of a non-central camera.
Upnp::CostParameters Upnp(const std::vector<Eigen::Vector3d>& ray_origins,
                          const std::vector<Eigen::Vector3d>& ray_directions,
                          const std::vector<Eigen::Vector3d>& world_points,
                          std::vector<Eigen::Quaterniond>* solution_rotations,
                          std::vector<Eigen::Vector3d>* solution_translations);

// Estimates the pose of a central camera.
Upnp::CostParameters Upnp(const std::vector<Eigen::Vector2d>& normalized_pixels,
                          const std::vector<Eigen::Vector3d>& world_points,
                          std::vector<Eigen::Quaterniond>* solution_rotations,
                          std::vector<Eigen::Vector3d>* solution_translations);

}  // namespace theia

#endif  // THEIA_SFM_POSE_UPNP_H_
