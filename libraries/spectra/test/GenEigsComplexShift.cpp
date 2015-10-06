#include <Eigen/Core>
#include <iostream>

#include <GenEigsSolver.h>
#include <MatOp/DenseGenComplexShiftSolve.h>

using namespace Spectra;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXcd ComplexMatrix;
typedef Eigen::VectorXcd ComplexVector;

template <int SelectionRule>
void run_test(const Matrix &mat, int k, int m, double sigmar, double sigmai)
{
    DenseGenComplexShiftSolve<double> op(mat);
    GenEigsComplexShiftSolver<double, SelectionRule, DenseGenComplexShiftSolve<double>> eigs(&op, k, m, sigmar, sigmai);
    eigs.init();
    int nconv = eigs.compute();
    int niter = eigs.num_iterations();
    int nops = eigs.num_operations();

    REQUIRE( nconv > 0 );

    ComplexVector evals = eigs.eigenvalues();
    ComplexMatrix evecs = eigs.eigenvectors();

    ComplexMatrix err = mat * evecs - evecs * evals.asDiagonal();

    INFO( "nconv = " << nconv );
    INFO( "niter = " << niter );
    INFO( "nops = " << nops );
    INFO( "||AU - UD||_inf = " << err.array().abs().maxCoeff() );
    INFO( "||AU - UD||_2 colwise =" << err.colwise().norm() );
    REQUIRE( err.array().abs().maxCoeff() == Approx(0.0) );
}


void run_test_sets(const Matrix &A, int k, int m, double sigmar, double sigmai)
{
    SECTION( "Largest Magnitude" )
    {
        run_test<LARGEST_MAGN>(A, k, m, sigmar, sigmai);
    }
    SECTION( "Largest Real Part" )
    {
        run_test<LARGEST_REAL>(A, k, m, sigmar, sigmai);
    }
    SECTION( "Largest Imaginary Part" )
    {
        run_test<LARGEST_IMAG>(A, k, m, sigmar, sigmai);
    }
    SECTION( "Smallest Magnitude" )
    {
        run_test<SMALLEST_MAGN>(A, k, m, sigmar, sigmai);
    }
    SECTION( "Smallest Real Part" )
    {
        run_test<SMALLEST_REAL>(A, k, m, sigmar, sigmai);
    }
    SECTION( "Smallest Imaginary Part" )
    {
        run_test<SMALLEST_IMAG>(A, k, m, sigmar, sigmai);
    }
}

TEST_CASE("Eigensolver of general real matrix [10x10]", "[eigs_gen]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(10, 10);
    int k = 3;
    int m = 6;
    double sigmar = 2.0;
    double sigmai = 1.0;

    run_test_sets(A, k, m, sigmar, sigmai);
}

TEST_CASE("Eigensolver of general real matrix [100x100]", "[eigs_gen]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(100, 100);
    int k = 10;
    int m = 20;
    double sigmar = 2.0;
    double sigmai = 1.0;

    run_test_sets(A, k, m, sigmar, sigmai);
}

TEST_CASE("Eigensolver of general real matrix [1000x1000]", "[eigs_gen]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(1000, 1000);
    int k = 20;
    int m = 50;
    double sigmar = 2.0;
    double sigmai = 1.0;

    run_test_sets(A, k, m, sigmar, sigmai);
}
