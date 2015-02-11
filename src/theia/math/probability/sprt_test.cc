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

#include <math.h>
#include <vector>
#include "gtest/gtest.h"

#include "theia/solvers/estimator.h"
#include "theia/math/probability/sequential_probability_ratio.h"
#include "theia/util/random.h"

using std::vector;

namespace theia {
namespace {

struct Point {
  double x;
  double y;
  Point() {}
  Point(double _x, double _y) : x(_x), y(_y) {}
};

// y = mx + b
struct Line {
  double m;
  double b;
  Line() {}
  Line(double _m, double _b) : m(_m), b(_b) {}
};
}  // namespace

class LineEstimator : public Estimator<Point, Line> {
 public:
  LineEstimator() {}
  ~LineEstimator() {}

  double SampleSize() const { return 2; }
  bool EstimateModel(const vector<Point>& data,
                     std::vector<Line>* models) const {
    Line model;
    model.m = (data[1].y - data[0].y) / (data[1].x - data[0].x);
    model.b = data[1].y - model.m * data[1].x;
    models->push_back(model);
    return true;
  }

  double Error(const Point& point, const Line& line) const {
    double a = -1.0 * line.m;
    double b = 1.0;
    double c = -1.0 * line.b;

    return fabs(a * point.x + b * point.y + c) / (sqrt(pow(a * a + b * b, 2)));
  }
};

// TODO(cmsweeney): Make this test a verification (i.e. is the value coming out
// accurate?) instead of just a sanity check.
TEST(SPRTTest, CalculateSPRTDecisionThreshold) {
  double sigma = 0.05;
  double epsilon = 0.1;
  double decision_threshold = CalculateSPRTDecisionThreshold(sigma, epsilon);
  VLOG(0) << "Decision threshold: " << decision_threshold;

  // Test with change of values for timing.
  decision_threshold = CalculateSPRTDecisionThreshold(sigma, epsilon, 200, 3);
  VLOG(0) << "Decision threshold: " << decision_threshold;
}

TEST(SPRTTest, SequentialProbabilityRatioTestPass) {
  // Create a set of points along y=x with a small random pertubation.
  vector<Point> input_points;
  for (int i = 0; i < 10000; ++i) {
    double noise_x = RandDouble(-1, 1);
    double noise_y = RandDouble(-1, 1);
    input_points.push_back(Point(i + noise_x, i + noise_y));
  }
  // Test for the correct line.
  Line fitting_line(1.0, 0.0);
  LineEstimator estimator;
  // Error threshold to consider points inliers vs outliers.
  double error_thresh = 0.5;
  // Estimate type 1 error.
  double sigma = 0.05;
  // Estimate of inlier ratio.
  double epsilon = 0.6;

  // Calculate the decision threshold.
  double decision_threshold = CalculateSPRTDecisionThreshold(sigma, epsilon);

  std::vector<double> residuals =
      estimator.Residuals(input_points, fitting_line);

  // Output parameters of SPRT.
  int num_tested_points;
  double observed_inlier_ratio;

  // Execute SPRT with a line we expect to fit the data.
  bool sprt_success = SequentialProbabilityRatioTest(
      residuals, error_thresh, sigma, epsilon, decision_threshold,
      &num_tested_points, &observed_inlier_ratio);
  EXPECT_TRUE(sprt_success);
}

TEST(SPRTTest, SequentialProbabilityRatioTestFail) {
  // Create a set of points along y=x with a small random pertubation.
  vector<Point> input_points;
  for (int i = 0; i < 10000; ++i) {
    double noise_x = RandDouble(-1, 1);
    double noise_y = RandDouble(-1, 1);
    input_points.push_back(Point(i + noise_x, i + noise_y));
  }
  // Test for the correct line.
  Line fitting_line(1.0, 0.0);
  LineEstimator estimator;
  // Error threshold to consider points inliers vs outliers.
  double error_thresh = 0.5;
  // Estimate type 1 error.
  double sigma = 0.05;
  // Estimate of inlier ratio.
  double epsilon = 0.6;

  // Calculate the decision threshold.
  double decision_threshold = CalculateSPRTDecisionThreshold(sigma, epsilon);

  // Output parameters of SPRT.
  int num_tested_points;
  double observed_inlier_ratio;

  // Execute SPRT with a few lines that do not fit the data.
  Line not_fitting_line(-1.0, 50);
  std::vector<double> residuals =
      estimator.Residuals(input_points, not_fitting_line);

  bool sprt_success = SequentialProbabilityRatioTest(
      residuals, error_thresh, sigma, epsilon, decision_threshold,
      &num_tested_points, &observed_inlier_ratio);
  EXPECT_FALSE(sprt_success);

  not_fitting_line = Line(1.0, 10);
  residuals = estimator.Residuals(input_points, not_fitting_line);
  sprt_success = SequentialProbabilityRatioTest(
      residuals, error_thresh, sigma, epsilon, decision_threshold,
      &num_tested_points, &observed_inlier_ratio);
  EXPECT_FALSE(sprt_success);

  not_fitting_line = Line(2.0, 0);
  residuals = estimator.Residuals(input_points, not_fitting_line);
  sprt_success = SequentialProbabilityRatioTest(
      residuals, error_thresh, sigma, epsilon, decision_threshold,
      &num_tested_points, &observed_inlier_ratio);
  EXPECT_FALSE(sprt_success);
}

}  // namespace theia
