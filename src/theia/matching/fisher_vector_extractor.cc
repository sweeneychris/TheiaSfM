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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

extern "C" {
#include <vl/fisher.h>
#include <vl/gmm.h>
}

#include <Eigen/Core>
#include <glog/logging.h>
#include <memory>
#include <vector>

#include "theia/matching/fisher_vector_extractor.h"

namespace theia {
namespace {
Eigen::MatrixXf ConvertVectorOfFeaturesToMatrix(
    const std::vector<Eigen::VectorXf>& features) {
  Eigen::MatrixXf feature_table(features[0].size(), features.size());
  for (int i = 0; i < feature_table.cols(); i++) {
    feature_table.col(i) = features[i];
  }
  return feature_table;
}
}  // namespace

// Wrapper for the VlFeat gmm.
class FisherVectorExtractor::GaussianMixtureModel {
 public:
  // NOTE: We need the unique_ptr to use vlfeat's delete function for the GMM.
  GaussianMixtureModel(const int num_clusters)
      : num_clusters_(num_clusters), gmm_(nullptr, vl_gmm_delete) {}

  int num_clusters() const { return num_clusters_; }

  // Computes the Gaussian mixture model based on the input training features.
  bool Compute(const Eigen::MatrixXf& data_points) {
    CHECK(!data_points.hasNaN());

    // Create the GMM instance.
    gmm_.reset(vl_gmm_new(VL_TYPE_FLOAT, data_points.rows(), num_clusters_));
    if (!gmm_) {
      return false;
    }
    // Compute the model.
    const double distortion =
        vl_gmm_cluster(gmm_.get(), data_points.data(), data_points.cols());
    CHECK(!std::isnan(distortion));
    return true;
  }

  // Helper accessor methods for retreiving VLFeat types.
  void const* GetMeans() const { return vl_gmm_get_means(gmm_.get()); }
  void const* GetCovariances() const {
    return vl_gmm_get_covariances(gmm_.get());
  }
  void const* GetPriors() const { return vl_gmm_get_priors(gmm_.get()); }

 private:
  // Number of clusters.
  const int num_clusters_;
  // The VlFeat GMM model.
  std::unique_ptr<VlGMM, void (*)(VlGMM*)> gmm_;
};

FisherVectorExtractor::FisherVectorExtractor(const Options& options)
    : gmm_(new GaussianMixtureModel(options.num_gmm_clusters)),
      training_feature_sampler_(options.max_num_features_for_training) {}

FisherVectorExtractor::~FisherVectorExtractor() {}

void FisherVectorExtractor::AddFeaturesForTraining(
    const std::vector<Eigen::VectorXf>& features) {
  for (const Eigen::VectorXf& feature : features) {
    CHECK(!feature.hasNaN()) << "Feature: " << feature.transpose();
    training_feature_sampler_.AddElementToSampler(feature);
  }
}

bool FisherVectorExtractor::Train() {
  // Get the features randomly sampled for training.
  const auto& sampled_features = training_feature_sampler_.GetAllSamples();
  CHECK_GT(sampled_features.size(), 0);
  LOG(INFO) << "Training GMM for Fisher Vector extractin with "
            << sampled_features.size() << " features sampled from "
            << training_feature_sampler_.NumElementsAdded()
            << " total features.";

  // Train the GMM using the training feaures.
  const Eigen::MatrixXf feature_table =
      ConvertVectorOfFeaturesToMatrix(sampled_features);
  return gmm_->Compute(feature_table);
}

Eigen::VectorXf FisherVectorExtractor::ExtractGlobalDescriptor(
    const std::vector<Eigen::VectorXf>& features) {
  // Ensure there are input features and they are not zero dimensions.
  CHECK_GT(features.size(), 0);
  CHECK_GT(features[0].size(), 0);

  // Convert the features into a continuous memory block. The matrix is of size
  // D x N where D is the number of descrip
  const Eigen::MatrixXf feature_table =
      ConvertVectorOfFeaturesToMatrix(features);
  // Compute the fisher vector encoding.
  Eigen::VectorXf fisher_vector(2 * feature_table.rows() *
                                gmm_->num_clusters());
  vl_fisher_encode(fisher_vector.data(),
                   VL_TYPE_FLOAT,
                   gmm_->GetMeans(),
                   feature_table.rows(),
                   gmm_->num_clusters(),
                   gmm_->GetCovariances(),
                   gmm_->GetPriors(),
                   feature_table.data(),
                   feature_table.cols(),
                   VL_FISHER_FLAG_IMPROVED);
  DCHECK(std::isfinite(fisher_vector.sum()));
  return fisher_vector;
}

}  // namespace theia
