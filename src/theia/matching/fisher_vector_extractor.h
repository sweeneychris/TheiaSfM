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

#ifndef THEIA_MATCHING_FISHER_VECTOR_EXTRACTOR_H_
#define THEIA_MATCHING_FISHER_VECTOR_EXTRACTOR_H_

#include <Eigen/Core>
#include <memory>
#include <vector>

#include "theia/matching/global_descriptor_extractor.h"
#include "theia/math/reservoir_sampler.h"

namespace theia {

// A Fisher Vector is an image representation obtained by pooling local image
// features. It is frequently used as a global image descriptor in visual
// classification. A Gaussian Mixture Model is fitted to the training data, and
// the Fisher Vector is computed as as the mean and covariance deviation vectors
// from the modes of the distribution of the GMM.
class FisherVectorExtractor : public GlobalDescriptorExtractor {
 public:
  struct Options {
    // The number of cluster to use for the Gaussian Mixture Model that power
    // the Fisher Kernel.
    int num_gmm_clusters = 16;

    // If more than this number of features are added to the
    // FisherVectorExtractor then we randomly sample
    // max_num_features_for_training using a memory efficient Reservoir sampler
    // to avoid holding all features in memory.
    int max_num_features_for_training = 100000;
  };

  // The number of clusters to use for the GMM.
  FisherVectorExtractor(const Options& options);

  ~FisherVectorExtractor();

  // Add features to the descriptor extractor for training. This method may be
  // called multiple times to add multiple sets of features (e.g., once per
  // image) to the global descriptor extractor for training.
  void AddFeaturesForTraining(
      const std::vector<Eigen::VectorXf>& features) override;

  // Train the global descriptor extracto with the given set of feature
  // descriptors added with AddFeaturesForTraining. It is assumed that all
  // descriptors have the same length.
  bool Train() override;

  // Compute a global image descriptor for the set of input features.
  Eigen::VectorXf ExtractGlobalDescriptor(
      const std::vector<Eigen::VectorXf>& features) override;

 private:
  // A Gaussian Mixture Model is used to compute the Fisher Kernel.
  class GaussianMixtureModel;
  std::unique_ptr<GaussianMixtureModel> gmm_;

  // The GMM is trained from a set of feature descriptors. A reservoir sampler
  // is used to randomly sample features from an unknown number of input
  // features for training.
  ReservoirSampler<Eigen::VectorXf> training_feature_sampler_;
};

}  // namespace theia
#endif  // THEIA_MATCHING_FISHER_VECTOR_EXTRACTOR_H_
