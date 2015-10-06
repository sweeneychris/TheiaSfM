// The code was adapted from Eigen/src/Eigenvaleus/EigenSolver.h
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2010,2012 Jitse Niesen <jitse@maths.leeds.ac.uk>
// Copyright (C) 2015 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef UPPER_HESSENBERG_EIGEN_H
#define UPPER_HESSENBERG_EIGEN_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <stdexcept>

namespace Spectra {


template <typename Scalar = double>
class UpperHessenbergEigen
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    typedef Eigen::Ref<Matrix> GenericMatrix;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

    typedef std::complex<Scalar> Complex;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> ComplexMatrix;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> ComplexVector;

    int dim_n;                             // Size of the matrix
    Eigen::RealSchur<Matrix> m_realSchur;  // Schur decomposition solver
    Matrix m_matT;                         // Schur T matrix
    Matrix m_eivec;                        // Storing eigenvectors
    ComplexVector m_eivalues;              // Eigenvalues

    bool computed;

    // Complex scalar division
    static Complex cdiv(const Scalar& xr, const Scalar& xi, const Scalar& yr, const Scalar& yi)
    {
        Scalar r, d;
        if(std::abs(yr) > std::abs(yi))
        {
            r = yi/yr;
            d = yr + r*yi;
            return Complex((xr + r*xi)/d, (xi - r*xr)/d);
        } else {
            r = yr/yi;
            d = yi + r*yr;
            return Complex((r*xr + xi)/d, (r*xi - xr)/d);
        }
    }

    void doComputeEigenvectors()
    {
        const int size = m_eivec.cols();
        const Scalar eps = std::numeric_limits<Scalar>::epsilon();

        // inefficient! this is already computed in RealSchur
        Scalar norm(0);
        for(int j = 0; j < size; ++j)
        {
            norm += m_matT.row(j).segment((std::max)(j-1, 0), size-(std::max)(j-1, 0)).cwiseAbs().sum();
        }

        // Backsubstitute to find vectors of upper triangular form
        if(norm == Scalar(0))
            return;

        for(int n = size - 1; n >= 0; n--)
        {
            Scalar p = m_eivalues.coeff(n).real();
            Scalar q = m_eivalues.coeff(n).imag();

            // Scalar vector
            if(q == Scalar(0))
            {
                Scalar lastr(0), lastw(0);
                int l = n;

                m_matT.coeffRef(n,n) = 1.0;
                for(int i = n-1; i >= 0; i--)
                {
                    Scalar w = m_matT.coeff(i,i) - p;
                    Scalar r = m_matT.row(i).segment(l,n-l+1).dot(m_matT.col(n).segment(l, n-l+1));

                    if(m_eivalues.coeff(i).imag() < 0.0)
                    {
                        lastw = w;
                        lastr = r;
                    } else {
                        l = i;
                        if(m_eivalues.coeff(i).imag() == 0.0)
                        {
                            if (w != 0.0)
                                m_matT.coeffRef(i,n) = -r / w;
                            else
                                m_matT.coeffRef(i,n) = -r / (eps * norm);
                        }
                        else // Solve real equations
                        {
                            Scalar x = m_matT.coeff(i,i+1);
                            Scalar y = m_matT.coeff(i+1,i);
                            Scalar denom = (m_eivalues.coeff(i).real() - p) * (m_eivalues.coeff(i).real() - p) + m_eivalues.coeff(i).imag() * m_eivalues.coeff(i).imag();
                            Scalar t = (x * lastr - lastw * r) / denom;
                            m_matT.coeffRef(i,n) = t;
                            if(std::abs(x) > std::abs(lastw))
                                m_matT.coeffRef(i+1,n) = (-r - w * t) / x;
                            else
                                m_matT.coeffRef(i+1,n) = (-lastr - y * t) / lastw;
                        }

                        // Overflow control
                        Scalar t = std::abs(m_matT.coeff(i,n));
                        if((eps * t) * t > Scalar(1))
                            m_matT.col(n).tail(size-i) /= t;
                    }
                }
            } else if(q < Scalar(0) && n > 0) {  // Complex vector
                Scalar lastra(0), lastsa(0), lastw(0);
                int l = n-1;

                // Last vector component imaginary so matrix is triangular
                if(std::abs(m_matT.coeff(n,n-1)) > std::abs(m_matT.coeff(n-1,n)))
                {
                    m_matT.coeffRef(n-1,n-1) = q / m_matT.coeff(n,n-1);
                    m_matT.coeffRef(n-1,n) = -(m_matT.coeff(n,n) - p) / m_matT.coeff(n,n-1);
                }
                else
                {
                    Complex cc = cdiv(0.0, -m_matT.coeff(n-1,n), m_matT.coeff(n-1,n-1)-p, q);
                    m_matT.coeffRef(n-1,n-1) = cc.real();
                    m_matT.coeffRef(n-1,n) = cc.imag();
                }
                m_matT.coeffRef(n,n-1) = 0.0;
                m_matT.coeffRef(n,n) = 1.0;
                for(int i = n-2; i >= 0; i--)
                {
                    Scalar ra = m_matT.row(i).segment(l, n-l+1).dot(m_matT.col(n-1).segment(l, n-l+1));
                    Scalar sa = m_matT.row(i).segment(l, n-l+1).dot(m_matT.col(n).segment(l, n-l+1));
                    Scalar w = m_matT.coeff(i,i) - p;

                    if(m_eivalues.coeff(i).imag() < 0.0)
                    {
                        lastw = w;
                        lastra = ra;
                        lastsa = sa;
                    }
                    else
                    {
                        l = i;
                        if(m_eivalues.coeff(i).imag() == Scalar(0))
                        {
                            Complex cc = cdiv(-ra,-sa,w,q);
                            m_matT.coeffRef(i,n-1) = cc.real();
                            m_matT.coeffRef(i,n) = cc.imag();
                        }
                        else
                        {
                            // Solve complex equations
                            Scalar x = m_matT.coeff(i,i+1);
                            Scalar y = m_matT.coeff(i+1,i);
                            Scalar vr = (m_eivalues.coeff(i).real() - p) * (m_eivalues.coeff(i).real() - p) + m_eivalues.coeff(i).imag() * m_eivalues.coeff(i).imag() - q * q;
                            Scalar vi = (m_eivalues.coeff(i).real() - p) * Scalar(2) * q;
                            if((vr == 0.0) && (vi == 0.0))
                                vr = eps * norm * (std::abs(w) + std::abs(q) + std::abs(x) + std::abs(y) + std::abs(lastw));

                            Complex cc = cdiv(x*lastra-lastw*ra+q*sa, x*lastsa-lastw*sa-q*ra, vr, vi);
                            m_matT.coeffRef(i,n-1) = cc.real();
                            m_matT.coeffRef(i,n) = cc.imag();
                            if(std::abs(x) > (std::abs(lastw) + std::abs(q)))
                            {
                                m_matT.coeffRef(i+1,n-1) = (-ra - w * m_matT.coeff(i,n-1) + q * m_matT.coeff(i,n)) / x;
                                m_matT.coeffRef(i+1,n) = (-sa - w * m_matT.coeff(i,n) - q * m_matT.coeff(i,n-1)) / x;
                            }
                            else
                            {
                                cc = cdiv(-lastra-y*m_matT.coeff(i,n-1), -lastsa-y*m_matT.coeff(i,n), lastw, q);
                                m_matT.coeffRef(i+1,n-1) = cc.real();
                                m_matT.coeffRef(i+1,n) = cc.imag();
                            }
                        }

                        // Overflow control
                        Scalar t = std::max(std::abs(m_matT.coeff(i,n-1)), std::abs(m_matT.coeff(i,n)));
                        if((eps * t) * t > Scalar(1))
                            m_matT.block(i, n-1, size-i, 2) /= t;

                    }
                }

                // We handled a pair of complex conjugate eigenvalues, so need to skip them both
                n--;
            }
        }

        // Back transformation to get eigenvectors of original matrix
        Vector m_tmp(size);
        for(int j = size-1; j >= 0; j--)
        {
            m_tmp.noalias() = m_eivec.leftCols(j+1) * m_matT.col(j).segment(0, j+1);
            m_eivec.col(j) = m_tmp;
        }
    }

public:

    UpperHessenbergEigen() :
        dim_n(0), computed(false)
    {}

    UpperHessenbergEigen(ConstGenericMatrix &mat) :
        dim_n(mat.rows()), computed(false)
    {
        compute(mat);
    }

    void compute(ConstGenericMatrix &mat)
    {
        if(mat.rows() != mat.cols())
            throw std::invalid_argument("UpperHessenbergEigen: matrix must be square");

        dim_n = mat.rows();

        // Reduce to real Schur form
        Matrix Q = Matrix::Identity(dim_n, dim_n);
        m_realSchur.computeFromHessenberg(mat, Q, true);
        m_matT = m_realSchur.matrixT();
        m_eivec = m_realSchur.matrixU();

        // Compute eigenvalues from matT
        m_eivalues.resize(dim_n);
        int i = 0;
        while(i < dim_n)
        {
            // Real eigenvalue
            if(i == dim_n - 1 || m_matT.coeff(i+1, i) == Scalar(0))
            {
                m_eivalues.coeffRef(i) = m_matT.coeff(i, i);
                ++i;
            }
            else  // Complex eigenvalues
            {
                Scalar p = Scalar(0.5) * (m_matT.coeff(i, i) - m_matT.coeff(i+1, i+1));
                Scalar z = std::sqrt(std::abs(p * p + m_matT.coeff(i+1, i) * m_matT.coeff(i, i+1)));
                m_eivalues.coeffRef(i)   = Complex(m_matT.coeff(i+1, i+1) + p, z);
                m_eivalues.coeffRef(i+1) = Complex(m_matT.coeff(i+1, i+1) + p, -z);
                i += 2;
            }
        }

        // Compute eigenvectors
        doComputeEigenvectors();

        computed = true;
    }

    ComplexVector eigenvalues()
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergEigen: need to call compute() first");

        return m_eivalues;
    }

    ComplexMatrix eigenvectors()
    {
        if(!computed)
            throw std::logic_error("UpperHessenbergEigen: need to call compute() first");

        int n = m_eivec.cols();
        const Scalar prec = std::pow(std::numeric_limits<Scalar>::epsilon(), Scalar(2.0) / 3);

        ComplexMatrix matV(n, n);
        for(int j = 0; j < n; ++j)
        {
            if(std::abs(m_eivalues.coeff(j).imag()) <= prec || j + 1 == n)
            {
                // we have a real eigen value
                matV.col(j) = m_eivec.col(j).template cast<Complex>();
                matV.col(j).normalize();
            } else {
                // we have a pair of complex eigen values
                for(int i = 0; i < n; ++i)
                {
                    matV.coeffRef(i,j)   = Complex(m_eivec.coeff(i,j),  m_eivec.coeff(i,j+1));
                    matV.coeffRef(i,j+1) = Complex(m_eivec.coeff(i,j), -m_eivec.coeff(i,j+1));
                }
                matV.col(j).normalize();
                matV.col(j+1).normalize();
                ++j;
            }
        }

        return matV;
    }
};


} // namespace Spectra

#endif // UPPER_HESSENBERG_EIGEN_H
