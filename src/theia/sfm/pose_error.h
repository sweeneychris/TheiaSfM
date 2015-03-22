// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_POSE_ERROR_H_
#define THEIA_SFM_POSE_ERROR_H_

#include <glog/logging.h>
#include <memory>
#include <string>

#include "theia/math/histogram.h"
#include "theia/util/stringprintf.h"

namespace theia {

// A utility class for compute errors of camera pose. Repeatedly call AddError
// until all errors have been added then call MeanMedianHistogramMessage to get
// a string containing information about the pose errors.
class PoseError {
 public:
  // This histogram bins for pose errors must be specified.
  PoseError(const std::vector<double>& rotation_histogram_bins,
            const std::vector<double>& position_histogram_bins) {
    rotation_histogram_.reset(new Histogram<double>(rotation_histogram_bins));
    position_histogram_.reset(new Histogram<double>(position_histogram_bins));
  }
  void AddError(const double rotation_error, const double position_error) {
    rotation_error_.emplace_back(rotation_error);
    position_error_.emplace_back(position_error);

    rotation_histogram_->Add(rotation_error);
    position_histogram_->Add(position_error);
  }

  std::string PrintMeanMedianHistogram() {
    std::string message = "";
    if (rotation_error_.size() == 0) {
      message += "There were no  poses that were common to the model.";
      return message;
    }

    // Print rotation errors.
    std::sort(rotation_error_.begin(), rotation_error_.end());
    const double mean_rotation_error =
        std::accumulate(rotation_error_.begin(), rotation_error_.end(), 0.0) /
        static_cast<double>(rotation_error_.size());
    const double median_rotation_error =
        rotation_error_[rotation_error_.size() / 2];
    const std::string rotation_msg = rotation_histogram_->PrintString();
    message += StringPrintf(
        "Rotation Error:\nMean = %lf \nMedian = %lf\nHistogram:\n%s",
        mean_rotation_error, median_rotation_error, rotation_msg.c_str());


    // Print position errors.
    std::sort(position_error_.begin(), position_error_.end());
    const double mean_position_error =
        std::accumulate(position_error_.begin(), position_error_.end(), 0.0) /
        static_cast<double>(position_error_.size());
    const double median_position_error =
        position_error_[position_error_.size() / 2];
    const std::string position_msg = position_histogram_->PrintString();
    message += StringPrintf(
        "\nPosition Error:\nMean = %lf \nMedian = %lf\nHistogram:\n%s",
        mean_position_error, median_position_error, position_msg.c_str());
    return message;
  }

 private:
  std::unique_ptr<Histogram<double> > rotation_histogram_;
  std::unique_ptr<Histogram<double> > position_histogram_;

  std::vector<double> rotation_error_;
  std::vector<double> position_error_;
};

}  // namespace theia

#endif  // THEIA_SFM_POSE_ERROR_H_
