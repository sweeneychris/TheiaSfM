#include <Eigen/Core>
#include <iostream>

#include <SymEigsSolver.h>
#include <MatOp/DenseGenMatProd.h>

using namespace Spectra;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

template <int SelectionRule>
void run_test(const Matrix &mat, int k, int m)
{
    // Eigen::SelfAdjointEigenSolver<MatrixXd> eig(mat);
    // std::cout << "all eigenvalues = \n" << eig.eigenvalues().transpose() << "\n";

    DenseGenMatProd<double> op(mat);
    SymEigsSolver<double, SelectionRule, DenseGenMatProd<double>> eigs(&op, k, m);
    eigs.init();
    int nconv = eigs.compute();
    int niter = eigs.num_iterations();
    int nops = eigs.num_operations();

    REQUIRE( nconv > 0 );

    Vector evals = eigs.eigenvalues();
    Matrix evecs = eigs.eigenvectors();

    // std::cout << "computed eigenvalues D = \n" << evals.transpose() << "\n";
    // std::cout << "computed eigenvectors U = \n" << evecs << "\n\n";
    Matrix err = mat * evecs - evecs * evals.asDiagonal();

    INFO( "nconv = " << nconv );
    INFO( "niter = " << niter );
    INFO( "nops = " << nops );
    INFO( "||AU - UD||_inf = " << err.array().abs().maxCoeff() );
    REQUIRE( err.array().abs().maxCoeff() == Approx(0.0) );
}

void run_test_sets(const Matrix &mat, int k, int m)
{
    SECTION( "Largest Magnitude" )
    {
        run_test<LARGEST_MAGN>(mat, k, m);
    }
    SECTION( "Largest Value" )
    {
        run_test<LARGEST_ALGE>(mat, k, m);
    }
    SECTION( "Smallest Magnitude" )
    {
        run_test<SMALLEST_MAGN>(mat, k, m);
    }
    SECTION( "Smallest Value" )
    {
        run_test<SMALLEST_ALGE>(mat, k, m);
    }
    SECTION( "Both Ends" )
    {
        run_test<BOTH_ENDS>(mat, k, m);
    }
}

TEST_CASE("Eigensolver of symmetric real matrix [10x10]", "[eigs_sym]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(10, 10);
    Matrix M = A + A.transpose();
    int k = 3;
    int m = 6;

    run_test_sets(M, k, m);
}

TEST_CASE("Eigensolver of symmetric real matrix [100x100]", "[eigs_sym]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(100, 100);
    Matrix M = A + A.transpose();
    int k = 10;
    int m = 20;

    run_test_sets(M, k, m);
}

TEST_CASE("Eigensolver of symmetric real matrix [1000x1000]", "[eigs_sym]")
{
    srand(123);

    Matrix A = Eigen::MatrixXd::Random(1000, 1000);
    Matrix M = A + A.transpose();
    int k = 20;
    int m = 50;

    run_test_sets(M, k, m);
}
