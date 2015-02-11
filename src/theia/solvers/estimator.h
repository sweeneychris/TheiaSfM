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

#ifndef THEIA_SOLVERS_ESTIMATOR_H_
#define THEIA_SOLVERS_ESTIMATOR_H_

#include <glog/logging.h>
#ifdef THEIA_USE_OPENMP
#include <omp.h>
#endif
#include <vector>

namespace theia {
// Templated class for estimating a model for RANSAC. This class is purely a
// virtual class and should be implemented for the specific task that RANSAC is
// being used for. Two methods must be implemented: EstimateModel and Error. All
// other methods are optional, but will likely enhance the quality of the RANSAC
// output.
//
// NOTE: RANSAC, ARRSAC, and other solvers work best if Datum and Model are
// lightweight classes or structs.

template <typename DatumType, typename ModelType> class Estimator {
 public:
  typedef DatumType Datum;
  typedef ModelType Model;

  Estimator() {}
  virtual ~Estimator() {}

  // Get the minimum number of samples needed to generate a model.
  virtual double SampleSize() const = 0;

  // Given a set of data points, estimate the model. Users should implement this
  // function appropriately for the task being solved. Returns true for
  // successful model estimation (and outputs model), false for failed
  // estimation. Typically, this is a minimal set, but it is not required to be.
  virtual bool EstimateModel(const std::vector<Datum>& data,
                             std::vector<Model>* model) const = 0;

  // Estimate a model from a non-minimal sampling of the data. E.g. for a line,
  // use SVD on a set of points instead of constructing a line from two points.
  // By default, this simply implements the minimal case.
  virtual bool EstimateModelNonminimal(const std::vector<Datum>& data,
                                       std::vector<Model>* model) const {
    return EstimateModel(data, model);
  }

  // Refine the model based on an updated subset of data, and a pre-computed
  // model. Can be optionally implemented.
  virtual bool RefineModel(const std::vector<Datum>& data, Model* model) const {
    return true;
  }

  // Given a model and a data point, calculate the error. Users should implement
  // this function appropriately for the task being solved.
  virtual double Error(const Datum& data, const Model& model) const = 0;

  // Compute the residuals of many data points. By default this is just a loop
  // that calls Error() on each data point, but this function can be useful if
  // the errors of multiple points may be estimated simultanesously (e.g.,
  // matrix multiplication to compute the reprojection error of many points at
  // once).
  virtual std::vector<double> Residuals(const std::vector<Datum>& data,
                                        const Model& model) const {
    std::vector<double> residuals(data.size());
#pragma omp parallel for
    for (int i = 0; i < data.size(); i++) {
      residuals[i] = Error(data[i], model);
    }
    return residuals;
  }

  // Returns the set inliers of the data set based on the error threshold
  // provided.
  std::vector<int> GetInliers(const std::vector<Datum>& data,
                              const Model& model,
                              double error_threshold) const {
    std::vector<int> inliers;
    inliers.reserve(data.size());
    for (int i = 0; i < data.size(); i++) {
      if (Error(data[i], model) < error_threshold) {
        inliers.push_back(i);
      }
    }
    return inliers;
  }

  // Enable a quick check to see if the model is valid. This can be a geometric
  // check or some other verification of the model structure.
  virtual bool ValidModel(const Model& model) const { return true; }
};

}  // namespace theia

#endif  // THEIA_SOLVERS_ESTIMATOR_H_
