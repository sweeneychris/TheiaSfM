#include <Eigen/Core>
#include <iostream>

#include <GenEigsSolver.h>
#include <MatOp/DenseGenMatProd.h>

using namespace Spectra;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXcd ComplexMatrix;
typedef Eigen::VectorXcd ComplexVector;

template <int SelectionRule>
void run_test(const Matrix &mat, int k, int m)
{
    DenseGenMatProd<double> op(mat);
    GenEigsSolver<double, SelectionRule, DenseGenMatProd<double>> eigs(&op, k, m);
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
    REQUIRE( err.array().abs().maxCoeff() == Approx(0.0) );
}

void run_test_sets(const Matrix &A, int k, int m)
{
    SECTION( "Largest Magnitude" )
    {
        run_test<LARGEST_MAGN>(A, k, m);
    }
    SECTION( "Largest Real Part" )
    {
        run_test<LARGEST_REAL>(A, k, m);
    }
    SECTION( "Largest Imaginary Part" )
    {
        run_test<LARGEST_IMAG>(A, k, m);
    }
    SECTION( "Smallest Magnitude" )
    {
        run_test<SMALLEST_MAGN>(A, k, m);
    }
    SECTION( "Smallest Real Part" )
    {
        run_test<SMALLEST_REAL>(A, k, m);
    }
    SECTION( "Smallest Imaginary Part" )
    {
        run_test<SMALLEST_IMAG>(A, k, m);
    }
}

TEST_CASE("Eigensolver of general real matrix [10x10]", "[eigs_gen]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(10, 10);
    int k = 3;
    int m = 6;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of general real matrix [100x100]", "[eigs_gen]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(100, 100);
    int k = 10;
    int m = 20;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of general real matrix [1000x1000]", "[eigs_gen]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(1000, 1000);
    int k = 20;
    int m = 50;

    run_test_sets(A, k, m);
}
