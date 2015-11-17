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

    // These two are for convenience in adapting the tridiagonal_qr_step() function
    typedef Scalar RealScalar;
    typedef int Index;

    int n;
    Vector main_diag;     // Main diagonal elements of the matrix
    Vector sub_diag;      // Sub-diagonal elements of the matrix
    Matrix evecs;         // To store eigenvectors

    bool computed;

    static bool is_much_smaller_than(const Scalar &x, const Scalar &y,
        const Scalar &prec = Eigen::NumTraits<Scalar>::dummy_precision())
    {
        return Eigen::numext::abs2(x) <= Eigen::numext::abs2(y) * prec * prec;
    }

    // Adapted from Eigen/src/Eigenvaleus/SelfAdjointEigenSolver.h
    static void tridiagonal_qr_step(RealScalar *diag,
                                    RealScalar *subdiag, Index start,
                                    Index end, Scalar *matrixQ,
                                    Index n)
    {
        RealScalar td = (diag[end-1] - diag[end]) * RealScalar(0.5);
        RealScalar e = subdiag[end-1];
        // Note that thanks to scaling, e^2 or td^2 cannot overflow, however they can still
        // underflow thus leading to inf/NaN values when using the following commented code:
        //   RealScalar e2 = numext::abs2(subdiag[end-1]);
        //   RealScalar mu = diag[end] - e2 / (td + (td>0 ? 1 : -1) * sqrt(td*td + e2));
        // This explain the following, somewhat more complicated, version:
        RealScalar mu = diag[end];
        if(td == 0)
            mu -= std::abs(e);
        else
        {
            RealScalar e2 = Eigen::numext::abs2(subdiag[end-1]);
            RealScalar h = Eigen::numext::hypot(td, e);
            if(e2==0)  mu -= (e / (td + (td>0 ? 1 : -1))) * (e / h);
            else       mu -= e2 / (td + (td>0 ? h : -h));
        }

        RealScalar x = diag[start] - mu;
        RealScalar z = subdiag[start];
        for(Index k = start; k < end; ++k)
        {
            Eigen::JacobiRotation<RealScalar> rot;
            rot.makeGivens(x, z);

            // do T = G' T G
            RealScalar sdk = rot.s() * diag[k] + rot.c() * subdiag[k];
            RealScalar dkp1 = rot.s() * subdiag[k] + rot.c() * diag[k + 1];

            diag[k] = rot.c() * (rot.c() * diag[k] - rot.s() * subdiag[k]) - rot.s() * (rot.c() * subdiag[k] - rot.s() * diag[k + 1]);
            diag[k + 1] = rot.s() * sdk + rot.c() * dkp1;
            subdiag[k] = rot.c() * sdk - rot.s() * dkp1;

            if(k > start)
                subdiag[k - 1] = rot.c() * subdiag[k - 1] - rot.s() * z;

            x = subdiag[k];

            if(k < end - 1)
            {
                z = -rot.s() * subdiag[k+1];
                subdiag[k + 1] = rot.c() * subdiag[k + 1];
            }

            // apply the givens rotation to the unit matrix Q = Q * G
            if(matrixQ)
            {
                Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> > q(matrixQ, n, n);
                q.applyOnTheRight(k, k + 1, rot);
            }
        }
    }

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
                if(is_much_smaller_than(std::abs(subd[i]), (std::abs(maind[i]) + std::abs(maind[i + 1]))))
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

            tridiagonal_qr_step(maind, subd, start, end, evecs.data(), n);
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
