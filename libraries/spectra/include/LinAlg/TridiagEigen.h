// The code was adapted from Eigen/src/Eigenvaleus/SelfAdjointEigenSolver.h
//
// Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>
// Copyright (C) 2015 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TRIDIAG_EIGEN_H
#define TRIDIAG_EIGEN_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <stdexcept>

namespace Spectra {


template <typename Scalar = double>
class TridiagEigen
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    typedef Eigen::Ref<Matrix> GenericMatrix;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

    int n;
    Vector main_diag;     // Main diagonal elements of the matrix
    Vector sub_diag;      // Sub-diagonal elements of the matrix
    Matrix evecs;         // To store eigenvectors

    bool computed;

public:
    TridiagEigen() :
        n(0), computed(false)
    {}

    TridiagEigen(ConstGenericMatrix &mat) :
        n(mat.rows()), computed(false)
    {
        compute(mat);
    }

    void compute(ConstGenericMatrix &mat)
    {
        if(mat.rows() != mat.cols())
            throw std::invalid_argument("TridiagEigen: matrix must be square");

        n = mat.rows();
        main_diag = mat.diagonal();
        sub_diag = mat.diagonal(-1);
        evecs.resize(n, n);
        evecs.setIdentity();

        int end = n - 1;
        int start = 0;
        int iter = 0; // total number of iterations
        int info = 0;

        Scalar *maind = main_diag.data();
        Scalar *subd = sub_diag.data();

        while(end > 0)
        {
            for(int i = start; i < end; i++)
                if(Eigen::internal::isMuchSmallerThan(std::abs(subd[i]), (std::abs(maind[i]) + std::abs(maind[i + 1]))))
                    subd[i] = 0;

            // find the largest unreduced block
            while(end > 0 && subd[end - 1] == 0)
                end--;

            if(end <= 0)
                break;

            // if we spent too many iterations, we give up
            iter++;
            if(iter > 30 * n)
            {
                info = 1;
                break;
            }

            start = end - 1;
            while(start > 0 && subd[start - 1] != 0)
                start--;

            Eigen::internal::tridiagonal_qr_step<Eigen::ColMajor>(maind, subd, start, end, evecs.data(), n);
        }

        if(info > 0)
            throw std::logic_error("TridiagEigen: failed to compute all the eigenvalues");

        computed = true;
    }

    Vector eigenvalues()
    {
        if(!computed)
            throw std::logic_error("TridiagEigen: need to call compute() first");

        // After calling compute(), main_diag will contain the eigenvalues.
        return main_diag;
    }

    Matrix eigenvectors()
    {
        if(!computed)
            throw std::logic_error("TridiagEigen: need to call compute() first");

        return evecs;
    }
};


} // namespace Spectra

#endif // TRIDIAG_EIGEN_H
