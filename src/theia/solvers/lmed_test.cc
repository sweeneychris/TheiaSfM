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
// Author: Victor Fragoso (victor.fragoso@mail.wvu.edu)

#include <algorithm>
#include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "theia/solvers/estimator.h"
#include "theia/solvers/lmed_quality_measurement.h"
#include "theia/util/random.h"

namespace theia {
namespace {
// Number of synthetic points. Note that the points have a 0.5 contamination
// rate, which is the maximum contamination rate that LMeds can take.
constexpr int kNumInlierPoints = 5000;
constexpr int kNumOutlierPoints = 5000;

// TODO(vfragoso): These classes below  (Point, Line, and LineEstimator) can be
// put in a single file. Several tests such as ransac_test.cc and prosac_test.cc
// use these classes.
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

class LineEstimator : public Estimator<Point, Line> {
 public:
  LineEstimator() {}
  ~LineEstimator() {}

  double SampleSize() const { return 2; }
  bool EstimateModel(const std::vector<Point>& data,
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
    return fabs(a * point.x + b * point.y + c) / sqrt(a * a + b * b);
  }
};

class LmedTest : public ::testing::Test {
 public:
  static void SetUpTestCase() {
    input_points = new std::vector<Point>;
    input_points->reserve(kNumInlierPoints + kNumOutlierPoints);
    for (int i = 0; i < kNumInlierPoints; ++i) {
      input_points->emplace_back(i + RandGaussian(0.0, 0.1),
                                 i + RandGaussian(0.0, 0.1));
    }
    for (int i = 0; i < kNumOutlierPoints; ++i) {
      input_points->emplace_back(RandDouble(0.0, 10000),
                                 RandDouble(0.0, 10000));
    }
    // Reshuffle.
    std::random_shuffle(input_points->begin(), input_points->end());
  }

  static void TearDownTestCase() {
    delete input_points;
  }

  // Synthetic points.
  static std::vector<Point>* input_points;
};

std::vector<Point>* LmedTest::input_points = nullptr;

}  // namespace

// Tests the computation of the squared residuals by using the correct line
// model.
TEST_F(LmedTest, ComputingQualityMeasureOfCorrectModel) {
  LineEstimator line_estimator;
  LmedQualityMeasurement lmed_quality_measurement(line_estimator.SampleSize());
  Line correct_line(1.0, 0.0);
  std::vector<double> residuals(input_points->size());
  for (int i = 0; i < residuals.size(); ++i) {
    residuals[i] = line_estimator.Error(input_points->at(i), correct_line);
  }
  EXPECT_LT(lmed_quality_measurement.ComputeCost(residuals), 0.5);
  EXPECT_NEAR(lmed_quality_measurement.GetInlierRatio(), 0.5, 0.1);
}

// TODO(vfragoso): Add tests for the LMed estimator.

}  // namespace theia
