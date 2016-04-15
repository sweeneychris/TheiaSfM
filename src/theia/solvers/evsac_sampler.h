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
// Author: Victor Fragoso (vfragoso@cs.ucsb.edu)

#ifndef THEIA_SOLVERS_EVSAC_SAMPLER_H_
#define THEIA_SOLVERS_EVSAC_SAMPLER_H_

#include <Eigen/Core>
#include <glog/logging.h>
#include <optimo/solvers/primal_dual_qp.h>
#include <statx/distributions/evd/gev.h>
#include <statx/distributions/gamma.h>
#include <statx/distributions/rayleigh.h>
#include <statx/utils/ecdf.h>

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <random>
#include <vector>

#include "theia/solvers/sampler.h"
#include "theia/util/util.h"

namespace theia {
// Fitting method used for fitting distributions.
//   MLE: Maximum likelihood estimation.
//   QUANTILE_NLS: Quantile non-linear least squares estimation.
enum FittingMethod {MLE = 0, QUANTILE_NLS = 1};

// EVSAC sampler implemented according to "EVSAC: Accelerating Hypotheses
// Generation by Modeling Matching Scores using Extreme Value Theory" by
// V. Fragoso, P. Sen, S. Rodriguez, M. Turk. from Proc. IEEE ICCV 2013.
// The idea of EVSAC is to model the statistics of the minimum distances
// computed when using the Nearest-Neighbor feature matcher. These minimums
// are the values used for taking the decisions and computing correspondences.
// The statistical model assumes that there are two processes generating these
// minimum distances. One that generates distances for correct correspondences,
// modeled by a gamma distribution, and another process generating distances for
// incorrect correspondences, modeled by the generalized extreme value (GEV)
// distribution. EVSAC tries to find the parameters for these two distributions
// as well as to estimate the mixing parameter, which happens to be an estimate
// of the inlier ratio.
template <class Datum> class EvsacSampler : public Sampler<Datum> {
 public:
  // Params:
  // min_num_samples:  The minimum number of samples to produce.
  // sorted_distances:  The matrix containing k L2 sorted distances. The matrix
  //   has num. of query features as rows and k columns.
  // predictor_threshold:  The threshold used to decide correct or incorrect
  //   matches/correspondences. The recommended value is 0.65.
  // fitting_method:  The fiting method to use, i.e.,  MLE or QUANTILE_NLS.
  EvsacSampler(
      const int min_num_samples,
      const Eigen::MatrixXd& sorted_distances,
      const double predictor_threshold,
      FittingMethod fitting_method)
      : Sampler<Datum>(min_num_samples), sorted_distances_(sorted_distances),
        predictor_threshold_(predictor_threshold),
        fitting_method_(fitting_method) {
    CHECK_GT(predictor_threshold_, 0.0);
    CHECK_LE(predictor_threshold_, 1.0);
  }

  ~EvsacSampler(void) {}

  // Helper function to calculate the number of distances required for a good
  // prediction of good correspondences using MRRayleigh predictor. Returns the
  // estimated size of the tail to use for prediction.
  //
  // Params:
  // num_reference_features:  The number of reference features.
  // percentile:  The percentile used to compute the size of the distribution
  //   tail from num_reference_features.
  static inline
  int CalculateKSmallestDistances(const int num_reference_features,
                                  const double percentile) {
    CHECK_GE(percentile, 0.0);
    CHECK_LE(percentile, 1.0);
    return static_cast<int>(num_reference_features * percentile);
  }

  // This function predicts if a correspondence is correct or incorrect using
  // the Meta-Recognition Rayleigh algorithm for distances. The function returns
  // true if a correspondence is predicted as correct and false otherwise.
  //
  // Params:
  // sorted_distances:  A container holding a sorted vector of the distances
  //   between a query and every reference descriptor corresponding to the 5% of
  //   the data.
  // predictor_threshold:  Confidence threshold to declare a correspondence as
  //   correct one. This confidence is within the interval of 0 and 1.
  static inline bool MRRayleigh(const Eigen::RowVectorXd& sorted_distances,
                                const double predictor_threshold) {
    CHECK_GT(predictor_threshold, 0.0);
    CHECK_LE(predictor_threshold, 1.0);
    const std::vector<double> tail(
        sorted_distances.data() + 1,
        sorted_distances.data() + sorted_distances.size());
    // Fit distribution tail!
    const double sigma = statx::distributions::raylfit(tail);
    // Calculate belief of correctness.
    const double confidence =
        1.0 - statx::distributions::raylcdf(sorted_distances[0], sigma);
    return confidence >= predictor_threshold;
  }

  // Parameters of the mixture of distributions (Gamma + GEV).
  struct MixtureModelParams {
    // Gamma Parameters (for correspondences modeled to be correct).
    double k;  // shape.
    double theta;  // scale.

    // GEV Parameters (for correspondences estimated to be incorrect).
    double xi;  // tail.
    double sigma;  // scale.
    double mu;  // location.

    // Inlier ratio (Mixture parameter).
    double inlier_ratio;  // Estimated inlier ratio.
  };

  // Calculates the mixture model parameters and the correctness probabilities
  // for every row in distances.
  // Params:
  // sorted_distances:  The matrix containing k L2 sorted distances. The
  //   matrix has num. of query features as rows and k columns.
  // predictor_threshold:  The threshold used to decide correct or
  //   incorrect matches/correspondences. The recommended value is 0.65.
  // fitting_method:  The fitting method MLE or QUANTILE_NLS (see statx doc).
  //   The recommended fitting method is the MLE estimation.
  // mixture_parameters:  The parameters of the pixture model.
  // probabilities:  The computed probabilities for every correspondence using
  //   the estimated mixture model.
  // sampling_weights:  The computed weights for non-uniform sampling. Samples
  //   with low probabilities of being correct are suppressed for sampling in a
  //   RANSAC scheme.
  static bool CalculateMixtureModel(
      const Eigen::MatrixXd& sorted_distances,
      const double predictor_threshold,
      const FittingMethod fitting_method,
      MixtureModelParams* mixture_parameters,
      std::vector<float>* probabilities,
      std::vector<float>* sampling_weights);

  // Implementing the Sample method.
  bool Sample(
      const std::vector<Datum>& data, std::vector<Datum>* subset) override;

  // Implementing the Initialize method.
  bool Initialize() override;

 protected:
  // L2 descriptor sorted sorted_distances.
  // rows: num of query features.
  // cols: k-th smallest distances.
  const Eigen::MatrixXd& sorted_distances_;
  // Threshold for predictor (MR-Rayleigh).
  const double predictor_threshold_;
  // Fitting method.
  FittingMethod fitting_method_;
  // Correspondence sampler following the computed probabilities.
  std::unique_ptr<std::discrete_distribution<>> correspondence_sampler_;
  // RNG
  std::mt19937 rng_;
  // Mixture Model Params.
  MixtureModelParams mixture_model_params_;

 private:
  // Extracts information from the sorted_distances to estimate the parameters
  // of the gamma and GEV distribution. Returns an estimated inlier ratio from
  // predictions.
  // Params:
  //   sorted_distances:  The matrix containing the NN distances.
  //   predictor_threshold: Threshold used for predicting correct or incorrect
  //     correspondences.
  //   predicted_correct_correspondences:  The binary vector indicating which
  //     correspondence are labeled as correct correspondence.
  //   smallest_distances:  The smallest distances in sorted_distances, this is
  //     the first column of sorted_distances matrix.
  //   negated_second_smallest_distances:  Negated second smallest distances,
  //     this is the negated second column of sorted_distances matrix. These are
  //     use to estimate the reversed generalized Pareto distribution (GEV)
  //     parameters.
  //   predicted_correct_correspondences_distances:  The distances of those
  //     correspondences that were predicted as correct.
  static double ExtractDataForFittingDistributions(
      const Eigen::MatrixXd& sorted_distances,
      const double predictor_threshold,
      std::vector<bool>* predictions,
      std::vector<double>* smallest_distances,
      std::vector<double>* negated_second_smallest_distances,
      std::vector<double>* distances_from_predicted_correspondences);

  // Estimates the two parameters of the Gamma distribution. Returns true upon
  // success, and false otherwise.
  // Params:
  //   predicted_correct_correspondences_distances:  The distances of those
  //     correspondences that were predicted as correct.
  //   mixture_model_parameters:  The mixture model parameters structures. The
  //     estimated gamma parameters are stored in this structure.
  static inline bool FitGamma(
      const std::vector<double>& predicted_correct_correspondences_distances,
      MixtureModelParams* mixture_model_parameters);

  // Estimates the three parameters of the generalized extreme value (GEV)
  // distribution. The function returns true upon success, and false otherwise.
  // Params:
  //   fitting_method:  The method for estimating the GEV parameters. This can
  //     maximum likelihood estimation, or non-linear least squares quantile
  //     based estimation.
  //   distances:  The random variables used for the estimation.
  static inline bool FitGEV(
      const FittingMethod fitting_method,
      const std::vector<double>& distances,
      MixtureModelParams* mixture_model_parameters);

  // This function calculates the inlier ratio by solving the following
  // least-squares (LS) problem:
  //
  // minimize_x 0.5 * norm(y - Ax)^2
  //
  // subject to
  //
  // [1 1] * x = 1
  // [0 0]' <= x <= [inlier_ratio_upper_bound 1]',
  //
  // where x = [inlier_ratio (1 - inlier_ratio)]'; which is our vector of
  // unknowns.
  // This LS problem can be transformed into a QP problem:
  //
  // minimize_x 0.5 * x' * A' * A * x - 2 * y' * A * x + y' * y
  //
  // subject to
  //
  // [1 1] * x = 1
  // [0 0]' <= x <= [inlier_ratio_upper_bound 1]'.
  //
  // This quadratic program finds the best inlier ratio such that the difference
  // between the empirical cdf of the smallest distances and the cdf obtained by
  // our mixture model (gamma + GEV) is minimal and that is less than the upper
  // bound estimated from the predictions.
  //
  // Params:
  //   inlier_ratio_upper_bound:  The upper bound for estimating the inlier
  //     ratio by solving the proposed QP problem.
  //   smallest_distances:  The smallest distances when matching.
  //   mixture_model_params:  The compute parameters for the gamma and GEV
  //     distributions.
  static bool EstimateInlierRatio(const double inlier_ratio_upper_bound,
                                  const std::vector<double>& smallest_distances,
                                  MixtureModelParams* mixture_model_params);

  // Computes the posterior and weights values for every correspondences. The
  // probabilities are used to generate a discrete distribution over the
  // correspondences and sample from them to generate hypotheses.
  // Params:
  //   num_correspondences:  The number of correspondences.
  //   mixture_model_params:  The compute mixture model parameters.
  //   smallest_distances:  The smallest distances in sorted_distances, this is
  //     the first column of sorted_distances matrix.
  //   predictions:  The correctness predictions for every correspndences.
  //   probabilities:  The computed probabilities for the correspondences used
  //     for sampling.
  //   sampling_weights:  The computed weights for non-uniform sampling.
  static void ComputePosteriorAndWeights(
      const int num_correspondences,
      const MixtureModelParams& mixture_model_params,
      const std::vector<double>& smallest_distances,
      const std::vector<bool>& predictions,
      std::vector<float>* probabilities,
      std::vector<float>* sampling_weights);

  typedef optimo::solvers::PrimalDualQP<
    double,
    2 /* Num. Unknowns */,
    4 /* Num. Inequalities */,
    1 /* Num. Equalities */ > PrimalDualQP;

  // Sets the inequality constraints.
  // The matrix Ain considers the case that inlier_ratio and (1 - inlier_ratio)
  // must be positive. Moreover, we set in this matrix the upperbound
  // constraints. The inequality constraints can be written in matrix form as
  // follows:
  // | -1  0 |                           |             0            |
  // |  0 -1 | |  inlier_ratio    |  <=  |             0            |
  // |  1  0 | | 1 - inlier_ratio |      | inlier_ratio_upper_bound |
  // |  0  1 |                           |             1            |
  //
  //    Ain           x              <=               bin.
  //
  // Params:
  //   inlier_ratio_upper_bound:  The upper bound for estimating the inlier
  //     ratio by solving the proposed QP problem.
  //   qp_params:  The QP params where the inequality constraints will be set.
  static inline void SetInequalityConstraints(
      const double inlier_ratio_upper_bound,  PrimalDualQP::Params* qp_params);

  // Sets the equality constraints.
  // The matrix Aeq considers the case that the vector of unknowns, x, must sum
  // up to 1:
  // | 1  1 | |   inlier_ratio   | = [ 1 ]
  //          | 1 - inlier_ratio |
  //   Aeq            x            =  beq.
  //
  // Params:
  //   qp_params:  The QP params where the equality constraints will be set.
  static inline
  void SetEqualityConstraints(PrimalDualQP::Params* qp_params);

  // Building the quadratic cost function for the QP:
  //
  // J(x) = x' * Q * x + d' * x,
  //
  // where
  // Q = A' * A, and d = A ' * empirical_cdf.
  // Matrix A is formed as follows:
  //
  // A(:, 0) = gammacdf(empirical_cdf_support; gamma_parameters)
  // A(:, 1) = 1 - gev(empirical_cdf_support; gev_parameters),
  //
  // where A(:, i) means the (i + 1)-th column, and the empirical_cdf_support is
  // the support of the empirical distribution function of the smallest
  // distances. Note that we set the second column to the reversed GEV; this is
  // done because we are considering minimums from the random process that
  // generates distances for incorrect correspondences.
  //
  // Params:
  //   smallest_distances:  The smallest distances when matching.
  //   qp_params:  The QP params where the equality constraints will be set.
  static void SetQPCostFunctionParams(
      const std::vector<double>& smallest_distances,
      const MixtureModelParams& mixture_model_params,
      PrimalDualQP::Params* qp_params);

  // Solves the QP problem to estimate the inlier ratio.
  // Params:
  //   inlier_ratio_upper_bound:  The inlier ratio upper bound estimated from
  //     predictions.
  //   qp_params:  The QP parameters that define the problem to solve.
  //   mixture_model_params:  Structure that holds the estimated inlier ratio.
  static inline
  bool SolveQPProblem(
      const double inlier_ratio_upper_bound,
      const PrimalDualQP::Params& qp_params,
      MixtureModelParams* mixture_model_params);

  DISALLOW_COPY_AND_ASSIGN(EvsacSampler);
};

// -------------------------- Implementation -------------------------------
template <class Datum>
bool EvsacSampler<Datum>::SolveQPProblem(
    const double inlier_ratio_upper_bound,
    const PrimalDualQP::Params& qp_params,
    MixtureModelParams* mixture_model_params) {
  Eigen::Vector2d unknowns_vector;
  // Setting initial point, i.e., guessed solution.
  unknowns_vector(0) = inlier_ratio_upper_bound / 2.0;
  unknowns_vector(1) = 1 - unknowns_vector(0);
  PrimalDualQP qp_solver;
  qp_solver.options.max_iter_ = 100;
  double min_value;
  const optimo::solvers::TERMINATION_TYPE termination_type =
      qp_solver(qp_params, &unknowns_vector, &min_value);
  VLOG(2) << "estimated inlier ratio=" << unknowns_vector(0)
          << " termination_type: " << termination_type;
  mixture_model_params->inlier_ratio = unknowns_vector(0);
  return termination_type == optimo::solvers::SOLVED;
}

template <class Datum>
void EvsacSampler<Datum>::SetQPCostFunctionParams(
    const std::vector<double>& smallest_distances,
    const MixtureModelParams& mixture_model_params,
    PrimalDualQP::Params* qp_params) {
  std::vector<double> empirical_cdf_vector;
  std::vector<double> empirical_cdf_support;
  statx::utils::ecdf(
      smallest_distances, &empirical_cdf_vector, &empirical_cdf_support);
  Eigen::Map<Eigen::VectorXd> empirical_cdf(
      &empirical_cdf_vector[0], empirical_cdf_vector.size());

  // Calculate matrix A.
  Eigen::MatrixXd A(empirical_cdf_support.size(), 2);

  for (int i = 0; i < empirical_cdf_support.size(); i++) {
    A(i, 0) = statx::distributions::gammacdf(
        empirical_cdf_support[i], mixture_model_params.k,
        mixture_model_params.theta);
    A(i, 1) = 1.0 - statx::distributions::evd::gevcdf(
        -empirical_cdf_support[i], mixture_model_params.mu,
        mixture_model_params.sigma, mixture_model_params.xi);
  }
  CHECK_NOTNULL(qp_params)->Q = A.transpose() * A;
  qp_params->d = -1.0 * A.transpose() * empirical_cdf;
}

template <class Datum>
void EvsacSampler<Datum>::SetEqualityConstraints(
    PrimalDualQP::Params* qp_params) {
  CHECK_NOTNULL(qp_params)->Aeq.setConstant(1.0);
  qp_params->beq(0, 0) = 1.0;
}

template <class Datum>
void EvsacSampler<Datum>::SetInequalityConstraints(
    const double inlier_ratio_upper_bound, PrimalDualQP::Params* qp_params) {
  CHECK_NOTNULL(qp_params)->Ain.setConstant(0.0);
  qp_params->Ain(0, 0) = -1.0;
  qp_params->Ain(1, 1) = -1.0;
  qp_params->Ain(2, 0) = 1.0;
  qp_params->Ain(3, 1) = 1.0;

  qp_params->bin.setConstant(0.0);
  qp_params->bin(2) = inlier_ratio_upper_bound;
  qp_params->bin(3) = 1.0;
}

template <class Datum>
bool EvsacSampler<Datum>::EstimateInlierRatio(
    const double inlier_ratio_upper_bound,
    const std::vector<double>& smallest_distances,
    MixtureModelParams* mixture_model_params) {
  PrimalDualQP::Params qp_params;

  // Set inequality constraints.
  SetInequalityConstraints(inlier_ratio_upper_bound, &qp_params);

  // Comppute equality constraints.
  SetEqualityConstraints(&qp_params);

  // Set the bulding QP params.
  SetQPCostFunctionParams(
      smallest_distances, *mixture_model_params, &qp_params);

  // Solve the QP.
  return SolveQPProblem(
      inlier_ratio_upper_bound, qp_params, mixture_model_params);
}

template <class Datum>
void EvsacSampler<Datum>::ComputePosteriorAndWeights(
    const int num_correspondences,
    const MixtureModelParams& mixture_model_params,
    const std::vector<double>& smallest_distances,
    const std::vector<bool>& predictions,
    std::vector<float>* probabilities,
    std::vector<float>* sampling_weights) {
  CHECK_NOTNULL(probabilities)->resize(num_correspondences);
  CHECK_NOTNULL(sampling_weights)->resize(num_correspondences);
  for (int i = 0; i < num_correspondences; i++) {
    // Calculate posterior.
    const double gam_val =
        mixture_model_params.inlier_ratio * statx::distributions::gammapdf(
            smallest_distances[i], mixture_model_params.k,
            mixture_model_params.theta);
    const double gev_val =
        (1.0 - mixture_model_params.inlier_ratio) *
        statx::distributions::evd::gevpdf(
            smallest_distances[i], mixture_model_params.mu,
            mixture_model_params.sigma, mixture_model_params.xi);
    const double posterior = gam_val / (gam_val + gev_val);
    // Removing those matches that are likely to be incorrect.
    (*probabilities)[i] = static_cast<float>(posterior);
    (*sampling_weights)[i] =
        predictions[i] ? static_cast<float>(posterior) : 0.0f;
  }
}

template <class Datum>
bool EvsacSampler<Datum>::FitGamma(
    const std::vector<double>& predicted_correct_correspondences_distances,
    MixtureModelParams* mixture_model_parameters) {
  const bool gam_success = statx::distributions::gammafit(
      predicted_correct_correspondences_distances,
      &CHECK_NOTNULL(mixture_model_parameters)->k,
      &mixture_model_parameters->theta);
  VLOG(2) << "Gamma distribution: k=" << mixture_model_parameters->k
          << " theta=" << mixture_model_parameters->theta
          << " flag: " << gam_success;
  return gam_success;
}

template <class Datum>
bool EvsacSampler<Datum>::FitGEV(
    const FittingMethod fitting_method,
    const std::vector<double>& negated_second_smallest_distances,
    MixtureModelParams* mixture_model_parameters) {
  const statx::distributions::evd::FitType fitting_type =
      (fitting_method == MLE) ?
      statx::distributions::evd::MLE : statx::distributions::evd::QUANTILE_NLS;
  const bool gev_success = gevfit(
      negated_second_smallest_distances,
      &CHECK_NOTNULL(mixture_model_parameters)->mu,
      &mixture_model_parameters->sigma, &mixture_model_parameters->xi,
      fitting_type);
  VLOG(2) << "GEV distribution: mu=" << mixture_model_parameters->mu
          << " sigma=" << mixture_model_parameters->sigma
          << " xi=" << mixture_model_parameters->xi << " flag: " << gev_success;
  return gev_success;
}

template <class Datum>
double EvsacSampler<Datum>::ExtractDataForFittingDistributions(
    const Eigen::MatrixXd& sorted_distances,
    const double predictor_threshold,
    std::vector<bool>* predicted_correct_correspondences,
    std::vector<double>* smallest_distances,
    std::vector<double>* negated_second_smallest_distances,
    std::vector<double>* predicted_correct_correspondences_distances) {
  CHECK_NOTNULL(predicted_correct_correspondences_distances)->reserve(
      sorted_distances.rows());
  CHECK_NOTNULL(predicted_correct_correspondences)->resize(
      sorted_distances.rows());
  CHECK_NOTNULL(smallest_distances)->resize(sorted_distances.rows());
  CHECK_NOTNULL(negated_second_smallest_distances)->resize(
      sorted_distances.rows());
  double estimated_inlier_ratio = 0.0;
  for (int i = 0; i < sorted_distances.rows(); ++i) {
    // Copying the first column from sorted distances to estimate the parameters
    // of the gamma distribution.
    (*smallest_distances)[i] = sorted_distances(i, 0);
    // Copying the second column from sorted distances as negated numbers to
    // estimate the parameters of the reversed GEV distribution.
    (*negated_second_smallest_distances)[i] = -sorted_distances(i, 1);
    // Saving the predictions.
    (*predicted_correct_correspondences)[i] =
        MRRayleigh(sorted_distances.row(i), predictor_threshold);
    if ((*predicted_correct_correspondences)[i]) {
      predicted_correct_correspondences_distances->push_back(
          sorted_distances(i, 0));
      estimated_inlier_ratio += 1.0;
    }
  }
  return estimated_inlier_ratio / sorted_distances.rows();
}

template <class Datum>
bool EvsacSampler<Datum>::CalculateMixtureModel(
    const Eigen::MatrixXd& sorted_distances,
    const double predictor_threshold,
    const FittingMethod fitting_method,
    MixtureModelParams* mixture_model_params,
    std::vector<float>* probabilities,
    std::vector<float>* sampling_weights) {
  CHECK_NOTNULL(mixture_model_params);
  CHECK_NOTNULL(probabilities);
  using std::vector;

  // A container indicating if a query point is predicted as correct or
  // incorrect.
  vector<bool> predicted_correct_correspondences;
  // The minimum distances, i.e., the first column of sorted_distances.
  vector<double> smallest_distances;
  // The second column of sorted_distances.
  vector<double> negated_second_smallest_distances;
  // Distances of correspondences that were predicted as correct, i.e., inliers.
  vector<double> predicted_inlier_distances;
  const double inlier_ratio_upper_bound = ExtractDataForFittingDistributions(
      sorted_distances, predictor_threshold, &predicted_correct_correspondences,
      &smallest_distances, &negated_second_smallest_distances,
      &predicted_inlier_distances);

  // Fit gamma distribution.
  if (!FitGamma(predicted_inlier_distances, mixture_model_params)) {
    VLOG(2) << "FitGamma failed.";
    return false;
  }

  // Fit GEV distribution.
  if (!FitGEV(fitting_method, negated_second_smallest_distances,
              mixture_model_params)) {
    VLOG(2) << "FitGEV failed.";
    return false;
  }

  // Estimate inlier ratio.
  if (!EstimateInlierRatio(inlier_ratio_upper_bound, smallest_distances,
                           mixture_model_params)) {
    VLOG(2) << "EstimateInlierRatio failed.";
    return false;
  }

  // Calculate posterior and final weights.
  ComputePosteriorAndWeights(
      sorted_distances.rows(), *mixture_model_params,
      smallest_distances, predicted_correct_correspondences, probabilities,
      sampling_weights);

  return true;
}

template <class Datum>
bool EvsacSampler<Datum>::Initialize() {
  CHECK_GT(this->sorted_distances_.rows(), 0);
  CHECK_GT(this->sorted_distances_.cols(), 0);

  // Initialize RNG.
  std::random_device rd;
  rng_.seed(rd());

  // Calculate Mixture model.
  std::vector<float> probabilities;
  std::vector<float> sampling_weights;
  if (!CalculateMixtureModel(
          this->sorted_distances_, this->predictor_threshold_,
          this->fitting_method_, &this->mixture_model_params_,
          &probabilities, &sampling_weights)) {
    return false;
  }

  // Initialize sampler.
#if defined(_MSC_VER) & _MSC_VER < 1900
  // MSVC has one of the constructors missing. See
  // http://stackoverflow.com/questions/21959404/initialising-stddiscrete-distribution-in-vs2013
  std::size_t i(0);
  correspondence_sampler_.reset(
      new std::discrete_distribution<int>(sampling_weights.size(),
                                          0.0,  // dummy!
                                          0.0,  // dummy!
                                          [&sampling_weights, &i](double) {
                                            auto w = sampling_weights[i];
                                            ++i;
                                            return w;
                                          }));
#else
  correspondence_sampler_.reset(new std::discrete_distribution<int>(
      sampling_weights.begin(), sampling_weights.end()));
#endif

  return true;
}

template <class Datum>
bool EvsacSampler<Datum>::Sample(const std::vector<Datum>& data,
                                 std::vector<Datum>* subset) {
  CHECK_EQ(data.size(), sorted_distances_.rows());
  CHECK_NOTNULL(subset)->resize(this->min_num_samples_);
  std::vector<int> random_numbers;
  for (int i = 0; i < this->min_num_samples_; i++) {
    int rand_number;
    // Generate a random number that has not already been used.
    while (std::find(random_numbers.begin(), random_numbers.end(),
                     (rand_number = (*correspondence_sampler_)(rng_)))
           != random_numbers.end()) {}

    random_numbers.push_back(rand_number);
    subset->at(i) = data[rand_number];
  }
  return true;
}
}  // namespace theia

#endif  // THEIA_SOLVERS_EVSAC_SAMPLER_H_
