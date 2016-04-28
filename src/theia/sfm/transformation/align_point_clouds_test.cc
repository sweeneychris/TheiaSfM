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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/sfm/transformation/align_point_clouds.h"

namespace theia {
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::RowMajor;
using Eigen::Vector3d;

namespace {
double kEpsilon = 1e-6;

void UmeyamaSimpleTest() {
  std::vector<Vector3d> left = {
      Vector3d(0.4, -3.105, 2.147),
      Vector3d(1.293, 7.1982, -.068),
      Vector3d(-5.34, 0.708, -3.69),
      Vector3d(-.345, 1.987, 0.936),
      Vector3d(0.93, 1.45, 1.079),
      Vector3d(-3.15, -4.73, 2.49),
      Vector3d(2.401, -2.03, -1.87),
      Vector3d(3.192, -.573, 0.1),
      Vector3d(-2.53, 3.07, -5.19)};

  const Matrix3d rotation_mat =
      Eigen::AngleAxisd(DegToRad(15.0), Vector3d(1.0, -2.7, 1.9).normalized())
          .toRotationMatrix();
  const Vector3d translation_vec(0, 2, 2);
  const double expected_scale = 1.5;

  // Transform the points.
  std::vector<Vector3d> right;
  for (int i = 0; i < left.size(); i++) {
    Vector3d transformed_point =
        expected_scale * rotation_mat * left[i] + translation_vec;
    right.emplace_back(transformed_point);
  }

  // Compute the similarity transformation.
  Matrix3d rotation;
  Vector3d translation;
  double scale;
  AlignPointCloudsUmeyama(left, right, &rotation, &translation, &scale);

  // Ensure the calculated transformation is the same as the one we set.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_LT(std::abs(rotation(i, j) - rotation_mat(i, j)), kEpsilon);
    }
    ASSERT_LT(std::abs(translation(i) - translation_vec(i)), kEpsilon);
  }
  ASSERT_LT(fabs(expected_scale - scale), kEpsilon);
}

/// Check with the default implementation and the expected results
void UmeyamaWithWeigthsEqualToOne() {
    
    // The algorithm is fast so we can make several iterations
    const size_t nbIterations = 50;
    for(size_t iteration = 0; iteration < nbIterations; ++iteration){
        size_t nb = rand() % 1000 + 4; // 4 pts required
        std::vector<Vector3d> left(nb);
        std::vector<double> weights(nb, 1.0);
        
        for(size_t i = 0; i < nb; ++i)
            left[i] = Eigen::Vector3d::Random();
        
        const Matrix3d rotation_mat =
        Eigen::AngleAxisd(DegToRad(rand() % 360), Eigen::Vector3d::Random().normalized())
        .toRotationMatrix();
        const Vector3d translation_vec = Eigen::Vector3d::Random();
        // We try to avoid a scale of zero
        const double expected_scale = std::max(0.001, std::abs(1.0 +  (rand() % 2 == 0 ? 1.0 : -1.0) * (rand() % 100)/10.0));
        
        // Transform the points.
        std::vector<Vector3d> right;
        for (size_t i = 0; i < left.size(); i++) {
            Vector3d transformed_point =
            expected_scale * rotation_mat * left[i] + translation_vec;
            right.emplace_back(transformed_point);
        }
        
        // Compute the similarity transformation.
        Matrix3d rotationRef;
        Vector3d translationRef;
        double scaleRef;
        AlignPointCloudsUmeyama(left, right, &rotationRef, &translationRef, &scaleRef);
        
        // Ensure the calculated transformation is the same as the one we set.
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                ASSERT_LT(std::abs(rotationRef(i, j) - rotation_mat(i, j)), kEpsilon);
            }
            ASSERT_LT(std::abs(translationRef(i) - translation_vec(i)), kEpsilon);
        }
        ASSERT_LT(fabs(scaleRef - expected_scale), kEpsilon);
        
        
        Matrix3d rotation;
        Vector3d translation;
        double scale;
        AlignPointCloudsUmeyamaWithWeights(left, right, weights, &rotation, &translation, &scale);
        
        // Check with the expected result
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                ASSERT_LT(std::abs(rotation(i, j) - rotation_mat(i, j)), kEpsilon);
            }
            ASSERT_LT(std::abs(translation(i) - translation_vec(i)), kEpsilon);
        }
        ASSERT_LT(fabs(scale - expected_scale), kEpsilon);
        //
        
        // Check with the reference implementation
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                ASSERT_LT(std::abs(rotation(i, j) - rotationRef(i, j)), kEpsilon);
            }
            ASSERT_LT(std::abs(translation(i) - translationRef(i)), kEpsilon);
        }
        ASSERT_LT(fabs(scale - scaleRef), kEpsilon);
        //
    }
}
    
// It is not easy to find a formula for the weights that always works
// We need to make some statistics to see if the function works in most cases
//
// This test adds some noise on a given fraction of the points
// We try to estimate the weights from the errors of the first pass (setting weights to 1)
// We can expect that the errors are where the noise was added
// From this first result we can set a weight for each point.
// The weight is given by exp(-distance * distance)
// Where distance is the distance between right[i] and s * R * left[i] + T estimated by a first pass
// To see if the new estimation from these weights is better
// we need to check that the new parameters are closer to the expected parameters
void UmeyamaWithWeigthsAndNoise() {
    
    // Make some statistics
    const size_t nbIterations = 1000; // 1000 is not a problem the algorith is fast
    
    // 15 % of noise
    const float noiseRatio = 0.15f;
    
    // Percentage to consider the test succeeds (is it enough ?)
    const float percentageToConsiderThatTheTestSucceeds = 95.0f;
    
    size_t succeeded = 0;
    for(size_t iteration = 0; iteration < nbIterations; ++iteration){
        size_t nb = rand() % 1000 + 4; // 4 pts required
        std::vector<Vector3d> left(nb);
        std::vector<double> weights(nb, 1.0);
        
        for(size_t i = 0; i < nb; ++i)
            left[i] = Eigen::Vector3d::Random();
        
        const Matrix3d rotation_mat =
        Eigen::AngleAxisd(DegToRad(rand() % 360), Eigen::Vector3d::Random().normalized())
        .toRotationMatrix();
        const Vector3d translation_vec = Eigen::Vector3d::Random();
        const double expected_scale = std::max(0.001, std::abs(1.0 +  (rand() % 2 == 0 ? 1.0 : -1.0) * (rand() % 100)/10.0));
        
        // Transform the points.
        std::vector<Vector3d> right;
        for (size_t i = 0; i < left.size(); i++) {
            Vector3d transformed_point =
            expected_scale * rotation_mat * left[i] + translation_vec;
            right.emplace_back(transformed_point);
        }
        
        // Add noise on scale, point and translation
        for(size_t i = 0, end = (size_t)(noiseRatio * nb); i < end; ++i){
            size_t k = rand() % nb;
            double noiseOnScale = std::abs(expected_scale  + std::abs((rand() % 2 == 0 ? 1.0 : -1.0) * (rand() % 100)/1000.0));
            right[k] = noiseOnScale * rotation_mat * (left[k] + Vector3d::Random()) + translation_vec + Vector3d::Random();
        }
        
        
        // We need to find some weights
        Matrix3d rotationRef1;
        Vector3d translationRef1;
        double scaleRef1;
        AlignPointCloudsUmeyamaWithWeights(left, right, weights, &rotationRef1, &translationRef1, &scaleRef1);
        for(size_t i = 0; i < nb; ++i){
            double dist = (right[i] - (scaleRef1 * rotationRef1 * left[i] + translationRef1)).squaredNorm();
            weights[i]= exp(-dist);
        }
        //
        
        Matrix3d rotationRef2;
        Vector3d translationRef2;
        double scaleRef2;
        AlignPointCloudsUmeyamaWithWeights(left, right, weights, &rotationRef2, &translationRef2, &scaleRef2);
        
        
        // Check if the parameters are closer to real parameters ?
        bool conditionOnScale = (std::abs(scaleRef2 - expected_scale) < std::abs(scaleRef1 - expected_scale));
        bool conditionOnTranslation = (translation_vec-translationRef2).norm()< (translation_vec - translationRef1).norm();
        bool conditionOnRotation = (rotation_mat - rotationRef2).norm() < (rotation_mat - rotationRef1).norm();
        
        if (conditionOnScale && conditionOnTranslation && conditionOnRotation)
            ++succeeded;
    }
    
    ASSERT_LE(percentageToConsiderThatTheTestSucceeds, 100.0 * succeeded/(0.0 + nbIterations));
}

}  // namespace

TEST(AlignPointCloudsUmeyama, SimpleTest) {
  UmeyamaSimpleTest();
}
    
TEST(AlignPointCloudsUmeyamaWithWeights, CompareWithDefaultImplementation) {
    UmeyamaWithWeigthsEqualToOne();
}

TEST(AlignPointCloudsUmeyamaWithWeights, WeightsAndNoise) {
    UmeyamaWithWeigthsAndNoise();
}

}  // namespace theia
