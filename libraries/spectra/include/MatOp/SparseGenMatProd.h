// Copyright (C) 2015 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SPARSE_GEN_MAT_PROD_H
#define SPARSE_GEN_MAT_PROD_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Spectra {


///
/// \ingroup MatOp
///
/// This class defines the matrix-vector multiplication operation on a
/// sparse real matrix \f$A\f$, i.e., calculating \f$y=Ax\f$ for any vector
/// \f$x\f$. It is mainly used in the GenEigsSolver and
/// SymEigsSolver eigen solvers.
///
template <typename Scalar>
class SparseGenMatProd
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::SparseMatrix<Scalar> SparseMatrix;

    const SparseMatrix &mat;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat_ An **Eigen** sparse matrix object, whose type is
    /// `Eigen::SparseMatrix<Scalar, ...>`.
    ///
    SparseGenMatProd(SparseMatrix &mat_) :
        mat(mat_)
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

#endif // SPARSE_GEN_MAT_PROD_H
