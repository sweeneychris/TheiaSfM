// Test ../include/LinAlg/UpperHessenbergQR.h and
//      ../include/LinAlg/DoubleShiftQR.h
#include <Eigen/Core>
#include <Eigen/QR>
#include <LinAlg/UpperHessenbergQR.h>
#include <LinAlg/DoubleShiftQR.h>

using namespace Spectra;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

template <typename Solver, typename MatrixType>
void run_test(MatrixType &H)
{
    Solver decomp(H);
    int n = H.rows();

    // Obtain Q matrix
    MatrixXd I = MatrixXd::Identity(n, n);
    MatrixXd Q = I;
    decomp.apply_QY(Q);

    // Test orthogonality
    MatrixXd QtQ = Q.transpose() * Q;
    INFO( "||Q'Q - I||_inf = " << (QtQ - I).cwiseAbs().maxCoeff() );
    REQUIRE( (QtQ - I).cwiseAbs().maxCoeff() == Approx(0.0) );

    MatrixXd QQt = Q * Q.transpose();
    INFO( "||QQ' - I||_inf = " << (QQt - I).cwiseAbs().maxCoeff() );
    REQUIRE( (QQt - I).cwiseAbs().maxCoeff() == Approx(0.0) );

    // Calculate R = Q'H
    MatrixXd R = decomp.matrix_R();
    MatrixXd Rlower = R.triangularView<Eigen::Lower>();
    Rlower.diagonal().setZero();
    INFO( "Whether R is upper triangular, error = " << Rlower.cwiseAbs().maxCoeff() );
    REQUIRE( Rlower.cwiseAbs().maxCoeff() == Approx(0.0) );

    // Compare H and QR
    INFO( "||H - QR||_inf = " << (H - Q * R).cwiseAbs().maxCoeff() );
    REQUIRE( (H - Q * R).cwiseAbs().maxCoeff() == Approx(0.0) );

    // Test RQ
    MatrixXd rq = R;
    decomp.apply_YQ(rq);
    INFO( "max error of RQ = " << (decomp.matrix_RQ() - rq).cwiseAbs().maxCoeff() );
    REQUIRE( (decomp.matrix_RQ() - rq).cwiseAbs().maxCoeff() == Approx(0.0) );

    // Test "apply" functions
    MatrixXd Y = MatrixXd::Random(n, n);

    MatrixXd QY = Y;
    decomp.apply_QY(QY);
    INFO( "max error of QY = " << (QY - Q * Y).cwiseAbs().maxCoeff() );
    REQUIRE( (QY - Q * Y).cwiseAbs().maxCoeff() == Approx(0.0) );

    MatrixXd YQ = Y;
    decomp.apply_YQ(YQ);
    INFO( "max error of YQ = " << (YQ - Y * Q).cwiseAbs().maxCoeff() );
    REQUIRE( (YQ - Y * Q).cwiseAbs().maxCoeff() == Approx(0.0) );

    MatrixXd QtY = Y;
    decomp.apply_QtY(QtY);
    INFO( "max error of Q'Y = " << (QtY - Q.transpose() * Y).cwiseAbs().maxCoeff() );
    REQUIRE( (QtY - Q.transpose() * Y).cwiseAbs().maxCoeff() == Approx(0.0) );

    MatrixXd YQt = Y;
    decomp.apply_YQt(YQt);
    INFO( "max error of YQ' = " << (YQt - Y * Q.transpose()).cwiseAbs().maxCoeff() );
    REQUIRE( (YQt - Y * Q.transpose()).cwiseAbs().maxCoeff() == Approx(0.0) );

    // Test "apply" functions for vectors
    VectorXd y = VectorXd::Random(n);

    VectorXd Qy = y;
    decomp.apply_QY(Qy);
    INFO( "max error of Qy = " << (Qy - Q * y).cwiseAbs().maxCoeff() );
    REQUIRE( (Qy - Q * y).cwiseAbs().maxCoeff() == Approx(0.0) );

    VectorXd Qty = y;
    decomp.apply_QtY(Qty);
    INFO( "max error of Q'y = " << (Qty - Q.transpose() * y).cwiseAbs().maxCoeff() );
    REQUIRE( (Qty - Q.transpose() * y).cwiseAbs().maxCoeff() == Approx(0.0) );
}

TEST_CASE("QR of upper Hessenberg matrix", "[QR]")
{
    srand(123);
    int n = 100;
    MatrixXd m = MatrixXd::Random(n, n);
    m.array() -= 0.5;
    MatrixXd H = m.triangularView<Eigen::Upper>();
    H.diagonal(-1) = m.diagonal(-1);

    run_test< UpperHessenbergQR<double> >(H);

    MapMat Hmap(H.data(), H.rows(), H.cols());
    run_test< UpperHessenbergQR<double> >(Hmap);
}

TEST_CASE("QR of Tridiagonal matrix", "[QR]")
{
    srand(123);
    int n = 100;
    MatrixXd m = MatrixXd::Random(n, n);
    m.array() -= 0.5;
    MatrixXd H = MatrixXd::Zero(n, n);
    H.diagonal() = m.diagonal();
    H.diagonal(-1) = m.diagonal(-1);
    H.diagonal(1) = m.diagonal(-1);

    run_test< TridiagQR<double> >(H);

    MapMat Hmap(H.data(), H.rows(), H.cols());
    run_test< TridiagQR<double> >(Hmap);
}

TEST_CASE("QR decomposition with double shifts", "QR")
{
    srand(123);
    int n = 100;
    MatrixXd m = MatrixXd::Random(n, n);
    m.array() -= 0.5;
    MatrixXd H = m.triangularView<Eigen::Upper>();
    H.diagonal(-1) = m.diagonal(-1);
    H(1, 0) = 0;  // Test for the case when sub-diagonal element is zero

    double s = 2, t = 3;

    MatrixXd M = H * H - s * H + t * MatrixXd::Identity(n, n);
    Eigen::HouseholderQR<MatrixXd> qr(M);
    MatrixXd Q0 = qr.householderQ();

    DoubleShiftQR<double> decomp(H, s, t);
    MatrixXd Q = MatrixXd::Identity(n, n);
    decomp.apply_YQ(Q);

    // Equal up to signs
    INFO( "max error of Q = " << (Q.cwiseAbs() - Q0.cwiseAbs()).cwiseAbs().maxCoeff() );
    REQUIRE( (Q.cwiseAbs() - Q0.cwiseAbs()).cwiseAbs().maxCoeff() == Approx(0.0) );

    // Test Q'HQ
    INFO( "max error of Q'HQ = " << (decomp.matrix_QtHQ() - Q.transpose() * H * Q).cwiseAbs().maxCoeff() );
    REQUIRE( (decomp.matrix_QtHQ() - Q.transpose() * H * Q).cwiseAbs().maxCoeff() == Approx(0.0) );

    // Test apply functions
    VectorXd y = VectorXd::Random(n);
    MatrixXd Y = MatrixXd::Random(n / 2, n);

    VectorXd Qty = y;
    decomp.apply_QtY(Qty);
    INFO( "max error of Q'y = " << (Qty - Q.transpose() * y).cwiseAbs().maxCoeff() );
    REQUIRE( (Qty - Q.transpose() * y).cwiseAbs().maxCoeff() == Approx(0.0) );

    MatrixXd YQ = Y;
    decomp.apply_YQ(YQ);
    INFO( "max error of YQ = " << (YQ- Y * Q).cwiseAbs().maxCoeff() );
    REQUIRE( (YQ- Y * Q).cwiseAbs().maxCoeff() == Approx(0.0) );
}
