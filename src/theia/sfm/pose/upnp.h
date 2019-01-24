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
#include "theia/util/util.h"

namespace theia {
// This class computes the pose of a central or non-central camera using the
// Universal  Perspective N-point method from "UPnP: An Optimal O(n) Solution to
// the Absolute Pose Problem with Universal Applicability" by Kneip et al. Upnp
// solves for the pose by minimizing the following cost function:
//
// J(R, t) = \sum_i || depth_i * ray_direction_ + ray_origin_i - R * p_i -t||^2,
//
// where R and t are the rotation and translations, respectively. Upnp re-writes
// the cost function above as a function that only depends on the rotation
// matrix, yielding a cost function of the form:
//
// J(R) = vec(R)' * quadratic_penalty_matrix * vec(R) +
//        linear_penalty_vector * vec(R) + gamma,
//
// where vec(R) is a 10-vector formed as a function of the corresponding
// quaternion of R.
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
  Upnp(const bool use_minimal_template) :
      use_minimal_template_(use_minimal_template) {}
  Upnp() : Upnp(false) {}
  ~Upnp() = default;

  // The Upnp cost function can be rewritten as follows:
  //
  // J(R) = vec(R)' * quadratic_penatly_mat * vec(R) +
  //        2 * linear_penalty_vector' * vec(R) + gamma,
  //
  // where vec(R) is the vector shown in Eq. 12 of the Upnp paper.
  struct CostParameters {
    CostParameters() {
      // The quadratic penalty matrix.
      quadratic_penalty_mat.setZero();
      // The linear penalty matrix.
      linear_penalty_vector.setZero();
      // The constant term.
      gamma = 0.0;
    }
    ~CostParameters() = default;

    Matrix10d quadratic_penalty_mat;
    Vector10d linear_penalty_vector;
    double gamma;
  };

  // Estimates pose of a central and non-central camera. The function returns
  // true when a pose is found, and false otherwise.
  //
  // Params:
  //   ray_origins:  The origins of each of the ray directions.
  //   ray_directions:  The unit-vector representing the direction from the
  //     center of a camera to a 3D point.
  //   world_points:  The 3D points. For the i-th world point there must be a
  //     corresponding ray origin and direction at the i-th entry in ray_origins
  //     and ray_directions.
  //   solution_rotations:  The computed and candidate quaternions or rotations.
  //   solution_translations:  The estimated and candidate translations.
  //   solution_costs:  The costs of the solutions.
  bool EstimatePose(const std::vector<Eigen::Vector3d>& ray_origins,
                    const std::vector<Eigen::Vector3d>& ray_directions,
                    const std::vector<Eigen::Vector3d>& world_points,
                    std::vector<Eigen::Quaterniond>* solution_rotations,
                    std::vector<Eigen::Vector3d>* solution_translations,
                    std::vector<double>* solution_costs);

  bool EstimatePose(const std::vector<Eigen::Vector3d>& ray_origins,
                    const std::vector<Eigen::Vector3d>& ray_directions,
                    const std::vector<Eigen::Vector3d>& world_points,
                    std::vector<Eigen::Quaterniond>* solution_rotations,
                    std::vector<Eigen::Vector3d>* solution_translations) {
    return EstimatePose(ray_origins, ray_directions, world_points,
                        solution_rotations, solution_translations, nullptr);
  }

  // Getters.
  const CostParameters& cost_params() const {
    return cost_params_;
  }

  // Evaluates the Upnp cost function given a rotation and the cost parameters.
  // The evaluated cost is the following:
  //
  // J(R) = vec(R)' * quadratic_penatly_mat * vec(R) +
  //        2 * linear_penalty_vector' * vec(R) + gamma.
  //
  // The function returns the cost given the rotation and parameters.
  //
  // Params:
  //   parameters:  The cost parameters of the function (i.e., quadratic penalty
  //     matrix, linear penalty vector, and gamma.
  //   rotation:  The quaternion representing the rotation.
  static double EvaluateCost(const CostParameters& parameters,
                             const Eigen::Quaterniond& rotation);

  // Computes the Upnp cost residual for a 3D point given a rotation and
  // translation. The evaluated residual is the following:
  //
  // Residual = || depth * ray_direction + ray_origin - R * p - t||^2
  //
  // where R is the rotation, t is the translation, and p is the 3D point.
  // The function returns the computed residual.
  //
  // Params:
  //   ray_origin:  The origin of the ray direction.
  //   ray_direction:  The unit-vector representing the direction from the
  //     center of a camera to a 3D point.
  //   world_point:  The 3D points
  //   rotation:  The quaternion representing the rotation.
  //   translation:  The translation.
  static double ComputeResidual(const Eigen::Vector3d& ray_origin,
                                const Eigen::Vector3d& ray_direction,
                                const Eigen::Vector3d& world_point,
                                const Eigen::Quaterniond& rotation,
                                const Eigen::Vector3d& translation);

 private:
  // Cost parameters.
  CostParameters cost_params_;

  // Template matrices.
  RowMajorMatrixXd minimal_sample_template_matrix_;
  RowMajorMatrixXd non_minimal_sample_template_matrix_;

  // Flag that controls which polynomial solver to use. If set true, the minimal
  // samples will use the SolveForRotationsFromMinimalSample(), which is slow;
  // and uses SolveForRotationsFromNonMinimalSample() when set true, which is
  // faster than the former.
  const bool use_minimal_template_;

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

  DISALLOW_COPY_AND_ASSIGN(Upnp);
};

// Estimates the pose of a non-central camera. The function returns the cost
// parameters of the problem given the input datum.
//
// Params:
//   ray_origins:  The origins of each of the ray directions.
//   ray_directions:  The unit-vector representing the direction from the
//     center of a camera to a 3D point.
//   world_points:  The 3D points. For the i-th world point there must be a
//     corresponding ray origin and direction at the i-th entry in ray_origins
//     and ray_directions.
//   solution_rotations:  The computed and candidate quaternions or rotations.
//   solution_translations:  The estimated and candidate translations.
Upnp::CostParameters Upnp(const std::vector<Eigen::Vector3d>& ray_origins,
                          const std::vector<Eigen::Vector3d>& ray_directions,
                          const std::vector<Eigen::Vector3d>& world_points,
                          std::vector<Eigen::Quaterniond>* solution_rotations,
                          std::vector<Eigen::Vector3d>* solution_translations);

// Estimates the pose of a central camera. The function returns the cost
// parameters of the problem given the input datum.
//
// Params:
//   normalized_pixels:  The 2D pixel positions using the normalized image
//    coordinates.
//   world_points:  The 3D points. For the i-th world point there must be a
//     corresponding ray origin and direction at the i-th entry in ray_origins
//     and ray_directions.
//   solution_rotations:  The computed and candidate quaternions or rotations.
//   solution_translations:  The estimated and candidate translations.
Upnp::CostParameters Upnp(const std::vector<Eigen::Vector2d>& normalized_pixels,
                          const std::vector<Eigen::Vector3d>& world_points,
                          std::vector<Eigen::Quaterniond>* solution_rotations,
                          std::vector<Eigen::Vector3d>* solution_translations);

}  // namespace theia

#endif  // THEIA_SFM_POSE_UPNP_H_
