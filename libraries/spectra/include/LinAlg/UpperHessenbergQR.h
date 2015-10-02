// Copyright (C) 2015 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef UPPER_HESSENBERG_QR_H
#define UPPER_HESSENBERG_QR_H

#include <Eigen/Core>
#include <cmath>      // std::sqrt
#include <algorithm>  // std::fill, std::copy
#include <limits>     // std::numeric_limits
#include <stdexcept>  // std::logic_error

namespace Spectra {


///
/// \defgroup LinearAlgebra Linear Algebra
///
/// A number of classes for linear algebra operations.

///
/// \ingroup LinearAlgebra
///
/// Perform the QR decomposition of an upper Hessenberg matrix.
///
/// \tparam Scalar The element type of the matrix.
/// Currently supported types are `float`, `double` and `long double`.
///
template <typename Scalar = double>
class UpperHessenbergQR
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, 2, 2> Matrix22;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic> RowVector;
    typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;

    typedef Eigen::Ref<Matrix> GenericMatrix;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

protected:
    int n;
    Matrix mat_T;
    // Gi = [ cos[i]  sin[i]]
    //      [-sin[i]  cos[i]]
    // Q = G1 * G2 * ... * G_{n-1}
    Array rot_cos;
    Array rot_sin;
    bool computed;
public:
    ///
    /// Default constructor. Computation can
    /// be performed later by calling the compute() method.
    ///
    UpperHessenbergQR() :
        n(0), computed(false)
    {}

    ///
    /// Constructor to create an object that performs and stores the
    /// QR decomposition of an upper Hessenberg matrix `mat`.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the upper triangular and the lower subdiagonal parts of
    /// the matrix are used.
    ///
    UpperHessenbergQR(ConstGenericMatrix &mat) :
        n(mat.rows()),
        mat_T(n, n),
        rot_cos(n - 1),
        rot_sin(n - 1),
        computed(false)
    {
        compute(mat);
    }

    ///
    /// Conduct the QR factorization of an upper Hessenberg matrix.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the upper triangular and the lower subdiagonal parts of
    /// the matrix are used.
    ///
    virtual void compute(ConstGenericMatrix &mat)
    {
        n = mat.rows();
        mat_T.resize(n, n);
        rot_cos.resize(n - 1);
        rot_sin.resize(n - 1);

        std::copy(mat.data(), mat.data() + mat.size(), mat_T.data());

        Scalar xi, xj, r, c, s, eps = std::numeric_limits<Scalar>::epsilon();
        Scalar *Tii, *ptr;
        for(int i = 0; i < n - 1; i++)
        {
            Tii = &mat_T(i, i);

            // Make sure mat_T is upper Hessenberg
            // Zero the elements below mat_T(i + 1, i)
            std::fill(Tii + 2, Tii + n - i, Scalar(0));

            xi = Tii[0];  // mat_T(i, i)
            xj = Tii[1];  // mat_T(i + 1, i)
            r = std::sqrt(xi * xi + xj * xj);
            if(r <= eps)
            {
                rot_cos[i] = c = 1;
                rot_sin[i] = s = 0;
            } else {
                rot_cos[i] = c = xi / r;
                rot_sin[i] = s = -xj / r;
            }
            // For a complete QR decomposition,
            // we first obtain the rotation matrix
            // G = [ cos  sin]
            //     [-sin  cos]
            // and then do T[i:(i + 1), i:(n - 1)] = G' * T[i:(i + 1), i:(n - 1)]

            // Gt << c, -s, s, c;
            // mat_T.block(i, i, 2, n - i) = Gt * mat_T.block(i, i, 2, n - i);
            Tii[0] = r;    // mat_T(i, i)     => r
            Tii[1] = 0;    // mat_T(i + 1, i) => 0
            ptr = Tii + n; // mat_T(i, k), k = i+1, i+2, ..., n-1
            for(int j = i + 1; j < n; j++, ptr += n)
            {
                Scalar tmp = ptr[0];
                ptr[0] = c * tmp - s * ptr[1];
                ptr[1] = s * tmp + c * ptr[1];
            }

            // If we do not need to calculate the R matrix, then
            // only the cos and sin sequences are required.
            // In such case we only update T[i + 1, (i + 1):(n - 1)]
            // mat_T.block(i + 1, i + 1, 1, n - i - 1) *= c;
            // mat_T.block(i + 1, i + 1, 1, n - i - 1) += s * mat_T.block(i, i + 1, 1, n - i - 1);
        }

        computed = true;
    }

    ///
    /// Return the \f$R\f$ matrix in the QR decomposition, which is an
    /// upper triangular matrix.
    ///
    /// \return Returned matrix type will be `Eigen::Matrix<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    Matrix matrix_R()
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        return mat_T;
    }

    ///
    /// Return the \f$RQ\f$ matrix, the multiplication of \f$R\f$ and \f$Q\f$,
    /// which is an upper Hessenberg matrix.
    ///
    /// \return Returned matrix type will be `Eigen::Matrix<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    virtual Matrix matrix_RQ()
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        // Make a copy of the R matrix
        Matrix RQ = mat_T.template triangularView<Eigen::Upper>();

        Scalar *c = rot_cos.data(),
               *s = rot_sin.data();
        for(int i = 0; i < n - 1; i++)
        {
            // RQ[, i:(i + 1)] = RQ[, i:(i + 1)] * Gi
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]

            Scalar *Yi, *Yi1;
            Yi = &RQ(0, i);
            Yi1 = Yi + n;  // RQ(0, i + 1)
            for(int j = 0; j < i + 2; j++)
            {
                Scalar tmp = Yi[j];
                Yi[j]  = (*c) * tmp - (*s) * Yi1[j];
                Yi1[j] = (*s) * tmp + (*c) * Yi1[j];
            }

            /* Vector Yi = RQ.block(0, i, i + 2, 1);
            RQ.block(0, i, i + 2, 1)     = (*c) * Yi - (*s) * RQ.block(0, i + 1, i + 2, 1);
            RQ.block(0, i + 1, i + 2, 1) = (*s) * Yi + (*c) * RQ.block(0, i + 1, i + 2, 1); */
            c++;
            s++;
        }

        return RQ;
    }

    ///
    /// Apply the \f$Q\f$ matrix to a vector \f$y\f$.
    ///
    /// \param Y A vector that will be overwritten by the matrix product \f$Qy\f$.
    ///
    /// Vector type can be `Eigen::Vector<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    // Y -> QY = G1 * G2 * ... * Y
    void apply_QY(Vector &Y)
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        Scalar tmp;
        for(int i = n - 2; i >= 0; i--)
        {
            // Y[i:(i + 1)] = Gi * Y[i:(i + 1)]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            tmp      = Y[i];
            Y[i]     =  rot_cos[i] * tmp + rot_sin[i] * Y[i + 1];
            Y[i + 1] = -rot_sin[i] * tmp + rot_cos[i] * Y[i + 1];
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to a vector \f$y\f$.
    ///
    /// \param Y A vector that will be overwritten by the matrix product \f$Q'y\f$.
    ///
    /// Vector type can be `Eigen::Vector<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    // Y -> Q'Y = G_{n-1}' * ... * G2' * G1' * Y
    void apply_QtY(Vector &Y)
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        Scalar tmp;
        for(int i = 0; i < n - 1; i++)
        {
            // Y[i:(i + 1)] = Gi' * Y[i:(i + 1)]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            tmp      = Y[i];
            Y[i]     = rot_cos[i] * tmp - rot_sin[i] * Y[i + 1];
            Y[i + 1] = rot_sin[i] * tmp + rot_cos[i] * Y[i + 1];
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$QY\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> QY = G1 * G2 * ... * Y
    void apply_QY(GenericMatrix Y)
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        Scalar *c = rot_cos.data() + n - 2,
               *s = rot_sin.data() + n - 2;
        RowVector Yi(Y.cols()), Yi1(Y.cols());
        for(int i = n - 2; i >= 0; i--)
        {
            // Y[i:(i + 1), ] = Gi * Y[i:(i + 1), ]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi  = Y.row(i);
            Yi1 = Y.row(i + 1);
            Y.row(i)     =  (*c) * Yi + (*s) * Yi1;
            Y.row(i + 1) = -(*s) * Yi + (*c) * Yi1;
            c--;
            s--;
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$Q'Y\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> Q'Y = G_{n-1}' * ... * G2' * G1' * Y
    void apply_QtY(GenericMatrix Y)
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        Scalar *c = rot_cos.data(),
               *s = rot_sin.data();
        RowVector Yi(Y.cols()), Yi1(Y.cols());
        for(int i = 0; i < n - 1; i++)
        {
            // Y[i:(i + 1), ] = Gi' * Y[i:(i + 1), ]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi = Y.row(i);
            Yi1 = Y.row(i + 1);
            Y.row(i)     = (*c) * Yi - (*s) * Yi1;
            Y.row(i + 1) = (*s) * Yi + (*c) * Yi1;
            c++;
            s++;
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$YQ\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> YQ = Y * G1 * G2 * ...
    void apply_YQ(GenericMatrix Y)
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        Scalar *c = rot_cos.data(),
               *s = rot_sin.data();
        /*Vector Yi(Y.rows());
        for(int i = 0; i < n - 1; i++)
        {
            // Y[, i:(i + 1)] = Y[, i:(i + 1)] * Gi
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi = Y.col(i);
            Y.col(i)     = (*c) * Yi - (*s) * Y.col(i + 1);
            Y.col(i + 1) = (*s) * Yi + (*c) * Y.col(i + 1);
            c++;
            s++;
        }*/
        Scalar *Y_col_i, *Y_col_i1;
        int nrow = Y.rows();
        for(int i = 0; i < n - 1; i++)
        {
            Y_col_i  = &Y(0, i);
            Y_col_i1 = &Y(0, i + 1);
            for(int j = 0; j < nrow; j++)
            {
                Scalar tmp = Y_col_i[j];
                Y_col_i[j]  = (*c) * tmp - (*s) * Y_col_i1[j];
                Y_col_i1[j] = (*s) * tmp + (*c) * Y_col_i1[j];
            }
            c++;
            s++;
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$YQ'\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> YQ' = Y * G_{n-1}' * ... * G2' * G1'
    void apply_YQt(GenericMatrix Y)
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        Scalar *c = rot_cos.data() + n - 2,
               *s = rot_sin.data() + n - 2;
        Vector Yi(Y.rows());
        for(int i = n - 2; i >= 0; i--)
        {
            // Y[, i:(i + 1)] = Y[, i:(i + 1)] * Gi'
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi = Y.col(i);
            Y.col(i)     =  (*c) * Yi + (*s) * Y.col(i + 1);
            Y.col(i + 1) = -(*s) * Yi + (*c) * Y.col(i + 1);
            c--;
            s--;
        }
    }
};



///
/// Perform the QR decomposition of a tridiagonal matrix, a special
/// case of upper Hessenberg matrices.
///
/// \tparam Scalar The element type of the matrix.
/// Currently supported types are `float`, `double` and `long double`.
///
template <typename Scalar = double>
class TridiagQR: public UpperHessenbergQR<Scalar>
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

public:
    ///
    /// Default constructor. Computation can
    /// be performed later by calling the compute() method.
    ///
    TridiagQR() :
        UpperHessenbergQR<Scalar>()
    {}

    ///
    /// Constructor to create an object that performs and stores the
    /// QR decomposition of a tridiagonal matrix `mat`.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the major- and sub- diagonal parts of
    /// the matrix are used.
    ///
    TridiagQR(ConstGenericMatrix &mat) :
        UpperHessenbergQR<Scalar>()
    {
        this->compute(mat);
    }

    ///
    /// Conduct the QR factorization of a tridiagonal matrix.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the major- and sub- diagonal parts of
    /// the matrix are used.
    ///
    void compute(ConstGenericMatrix &mat)
    {
        this->n = mat.rows();
        this->mat_T.resize(this->n, this->n);
        this->rot_cos.resize(this->n - 1);
        this->rot_sin.resize(this->n - 1);

        this->mat_T.setZero();
        this->mat_T.diagonal() = mat.diagonal();
        this->mat_T.diagonal(1) = mat.diagonal(-1);
        this->mat_T.diagonal(-1) = mat.diagonal(-1);

        // A number of pointers to avoid repeated address calculation
        Scalar *Tii = this->mat_T.data(),  // pointer to T[i, i]
               *ptr,                       // some location relative to Tii
               *c = this->rot_cos.data(),  // pointer to the cosine vector
               *s = this->rot_sin.data(),  // pointer to the sine vector
               r, tmp,
               eps = std::numeric_limits<Scalar>::epsilon();
        for(int i = 0; i < this->n - 2; i++)
        {
            // Tii[0] == T[i, i]
            // Tii[1] == T[i + 1, i]
            r = std::sqrt(Tii[0] * Tii[0] + Tii[1] * Tii[1]);
            if(r <= eps)
            {
                *c = 1;
                *s = 0;
            } else {
                *c =  Tii[0] / r;
                *s = -Tii[1] / r;
            }

            // For a complete QR decomposition,
            // we first obtain the rotation matrix
            // G = [ cos  sin]
            //     [-sin  cos]
            // and then do T[i:(i + 1), i:(i + 2)] = G' * T[i:(i + 1), i:(i + 2)]

            // Update T[i, i] and T[i + 1, i]
            // The updated value of T[i, i] is known to be r
            // The updated value of T[i + 1, i] is known to be 0
            Tii[0] = r;
            Tii[1] = 0;
            // Update T[i, i + 1] and T[i + 1, i + 1]
            // ptr[0] == T[i, i + 1]
            // ptr[1] == T[i + 1, i + 1]
            ptr = Tii + this->n;
            tmp = *ptr;
            ptr[0] = (*c) * tmp - (*s) * ptr[1];
            ptr[1] = (*s) * tmp + (*c) * ptr[1];
            // Update T[i, i + 2] and T[i + 1, i + 2]
            // ptr[0] == T[i, i + 2] == 0
            // ptr[1] == T[i + 1, i + 2]
            ptr += this->n;
            ptr[0] = -(*s) * ptr[1];
            ptr[1] *= (*c);

            // Move from T[i, i] to T[i + 1, i + 1]
            Tii += this->n + 1;
            // Increase c and s by 1
            c++;
            s++;


            // If we do not need to calculate the R matrix, then
            // only the cos and sin sequences are required.
            // In such case we only update T[i + 1, (i + 1):(i + 2)]
            // this->mat_T(i + 1, i + 1) = (*c) * this->mat_T(i + 1, i + 1) + (*s) * this->mat_T(i, i + 1);
            // this->mat_T(i + 1, i + 2) *= (*c);
        }
        // For i = n - 2
        r = std::sqrt(Tii[0] * Tii[0] + Tii[1] * Tii[1]);
        if(r <= eps)
        {
            *c = 1;
            *s = 0;
        } else {
            *c =  Tii[0] / r;
            *s = -Tii[1] / r;
        }
        Tii[0] = r;
        Tii[1] = 0;
        ptr = Tii + this->n;  // points to T[i - 2, i - 1]
        tmp = *ptr;
        ptr[0] = (*c) * tmp - (*s) * ptr[1];
        ptr[1] = (*s) * tmp + (*c) * ptr[1];

        this->computed = true;
    }

    ///
    /// Return the \f$RQ\f$ matrix, the multiplication of \f$R\f$ and \f$Q\f$,
    /// which is a tridiagonal matrix.
    ///
    /// \return Returned matrix type will be `Eigen::Matrix<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    Matrix matrix_RQ()
    {
        if(!this->computed)
            throw std::logic_error("TridiagQR: need to call compute() first");

        // Make a copy of the R matrix
        Matrix RQ(this->n, this->n);
        RQ.setZero();
        RQ.diagonal() = this->mat_T.diagonal();
        RQ.diagonal(1) = this->mat_T.diagonal(1);

        // [m11  m12] will point to RQ[i:(i+1), i:(i+1)]
        // [m21  m22]
        Scalar *m11 = RQ.data(), *m12, *m21, *m22,
               *c = this->rot_cos.data(),
               *s = this->rot_sin.data(),
               tmp;
        for(int i = 0; i < this->n - 1; i++)
        {
            m21 = m11 + 1;
            m12 = m11 + this->n;
            m22 = m12 + 1;
            tmp = *m21;

            // Update diagonal and the below-subdiagonal
            *m11 = (*c) * (*m11) - (*s) * (*m12);
            *m21 = (*c) * tmp - (*s) * (*m22);
            *m22 = (*s) * tmp + (*c) * (*m22);

            // Move m11 to RQ[i+1, i+1]
            m11  = m22;
            c++;
            s++;
        }

        // Copy the below-subdiagonal to above-subdiagonal
        RQ.diagonal(1) = RQ.diagonal(-1);

        return RQ;
    }
};


} // namespace Spectra

#endif // UPPER_HESSENBERG_QR_H
