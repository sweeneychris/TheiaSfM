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
#include <random>

#include "gtest/gtest.h"

#include "theia/solvers/estimator.h"
#include "theia/solvers/prosac.h"
#include "theia/test/test_utils.h"

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
}  // namespace

TEST(ProsacTest, LineFitting) {
  // Create a set of points along y=x with a small random pertubation.
  std::default_random_engine generator(90);
  std::normal_distribution<double> gauss_distribution(0.0, 0.5);
  std::normal_distribution<double> small_distribution(0.0, 0.05);
  const int num_points = 10000;
  std::vector<Point> input_points(num_points);
  std::vector<double> confidence(num_points);

  for (int i = 0; i < num_points; ++i) {
    if (i < 300) {
      double noise_x = small_distribution(generator);
      double noise_y = small_distribution(generator);
      input_points[i] = Point(i + noise_x, i + noise_y);
      confidence[i] = 0.95;
    } else {
      double noise_x = gauss_distribution(generator);
      double noise_y = gauss_distribution(generator);
      input_points[i] = Point(i + noise_x, i + noise_y);
      confidence[i] = 0.1;
    }
  }

  LineEstimator line_estimator;
  Line line;
  RansacParameters params;
  params.rng = std::make_shared<RandomNumberGenerator>(113);
  params.error_thresh = 0.5;
  Prosac<LineEstimator> prosac_line(params, line_estimator);
  prosac_line.Initialize();
  RansacSummary summary;
  prosac_line.Estimate(input_points, &line, &summary);
  ASSERT_LT(fabs(line.m - 1.0), 0.1);
}

}  // namespace theia
