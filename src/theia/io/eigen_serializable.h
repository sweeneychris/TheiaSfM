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
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#ifndef THEIA_IO_EIGEN_SERIALIZABLE_H_
#define THEIA_IO_EIGEN_SERIALIZABLE_H_

#include <cereal/cereal.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <Eigen/Dense>

namespace cereal {

// This code is really ugly, but it allows for Eigen types to be read from or
// written to disk using the Cereal IO library.
template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options,
          int _MaxRows, int _MaxCols>
inline typename std::enable_if<
    traits::is_output_serializable<BinaryData<_Scalar>, Archive>::value,
    void>::type
save(Archive& ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows,
                                _MaxCols> const& m) {
  int32_t rows = m.rows();
  int32_t cols = m.cols();
  ar(rows);
  ar(cols);
  ar(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options,
          int _MaxRows, int _MaxCols>
inline typename std::enable_if<
    traits::is_input_serializable<BinaryData<_Scalar>, Archive>::value,
    void>::type
load(Archive& ar,
     Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m) {
  int32_t rows;
  int32_t cols;
  ar(rows);
  ar(cols);
  m.resize(rows, cols);
  ar(binary_data(m.data(),
                 static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options,
          int _MaxRows, int _MaxCols>
inline typename std::enable_if<
    traits::is_output_serializable<BinaryData<_Scalar>, Archive>::value,
    void>::type
save(Archive& ar, Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows,
                                _MaxCols> const& m) {
  int32_t rows = m.rows();
  int32_t cols = m.cols();
  ar(rows);
  ar(cols);
  ar(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options,
          int _MaxRows, int _MaxCols>
inline typename std::enable_if<
    traits::is_input_serializable<BinaryData<_Scalar>, Archive>::value,
    void>::type
load(Archive& ar,
     Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m) {
  int32_t rows;
  int32_t cols;
  ar(rows);
  ar(cols);
  m.resize(rows, cols);
  ar(binary_data(m.data(),
                 static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
}

}  // namespace cereal

#endif  // THEIA_IO_EIGEN_SERIALIZABLE_H_
