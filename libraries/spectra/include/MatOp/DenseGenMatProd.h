// Copyright (C) 2015 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DENSE_GEN_MAT_PROD_H
#define DENSE_GEN_MAT_PROD_H

#include <Eigen/Core>

namespace Spectra {


///
/// \defgroup MatOp Matrix Operations
///
/// Define matrix operations on existing matrix objects
///

///
/// \ingroup MatOp
///
/// This class defines the matrix-vector multiplication operation on a
/// general real matrix \f$A\f$, i.e., calculating \f$y=Ax\f$ for any vector
/// \f$x\f$. It is mainly used in the GenEigsSolver and
/// SymEigsSolver eigen solvers.
///
template <typename Scalar>
class DenseGenMatProd
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Matrix> MapMat;
    typedef Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > MapVec;

    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

    const MapMat mat;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat_ An **Eigen** matrix object, whose type can be
    /// `Eigen::Matrix<Scalar, ...>` (e.g. `Eigen::MatrixXd` and
    /// `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    DenseGenMatProd(ConstGenericMatrix &mat_) :
        mat(mat_.data(), mat_.rows(), mat_.cols())
    {}

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    int rows() { return mat.rows(); }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    int cols() { return mat.cols(); }

    ///
    /// Perform the matrix-vector multiplication operation \f$y=Ax\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = A * x_in
    void perform_op(Scalar *x_in, Scalar *y_out)
    {
        MapVec x(x_in, mat.cols());
        MapVec y(y_out, mat.rows());
        y.noalias() = mat * x;
    }
};


} // namespace Spectra

#endif // DENSE_GEN_MAT_PROD_H
