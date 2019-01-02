// This file was created by Steffen Urban (urbste@googlemail.com) or company address (steffen.urban@zeiss.com)
// December 2018


#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <random>

#include "theia/util/random.h"
#include "theia/sfm/pose/four_point_focal_length_radial_distortion.h"

namespace theia
{

theia::RandomNumberGenerator rng(42);

double sgn(double val) {
  return (0.0 < val) - (val < 0.0);
}

void Solve(Eigen::Matrix<double, 64, 1> data,
    Eigen::Matrix<std::complex<double>, 5, 13>& sols)
{
    const double t2 = data(27) * data(27);
    const double t3 = data(29) * data(29);
    const double t4 = data(30) * data(30);
    const double t5 = data(19) * data(19);
    const double t6 = data(21) * data(21);
    const double t7 = data(22) * data(22);
    const double t8 = data(11) * data(11);
    const double t9 = data(13) * data(13);
    const double t10 = data(14) * data(14);
    const double t11 = data(3) * data(3);
    const double t12 = data(5) * data(5);
    const double t13 = data(6) * data(6);
    const double t14 = data(26) * data(26);
    const double t15 = data(31) * data(31);
    const double t16 = data(18)*data(26)*2.0;
    const double t17 = data(18) * data(18);
    const double t18 = data(23) * data(23);
    const double t19 = data(10)*data(26)*2.0;
    const double t20 = data(10)*data(18)*2.0;
    const double t21 = data(10) * data(10);
    const double t22 = data(15) * data(15);
    const double t23 = data(2)*data(26)*2.0;
    const double t24 = data(2)*data(18)*2.0;
    const double t25 = data(2)*data(10)*2.0;
    const double t26 = data(2) * data(2);
    const double t27 = data(7) * data(7);
    Eigen::Matrix<double, 328, 2> coeffs;
    Eigen::Matrix<double, 327, 1> coeffs1;
    coeffs(1, 1) = data(29)*data(57) + data(30)*data(58);
    coeffs(2, 1) = data(31) + data(29)*data(54) + data(30)*data(55);
    coeffs(3, 1) = data(29)*data(51) + data(30)*data(52);
    coeffs(4, 1) = data(29)*data(39) + data(30)*data(40) + data(21)*data(57) + data(22)*data(58);
    coeffs(5, 1) = data(23) + data(21)*data(54) + data(22)*data(55);
    coeffs(6, 1) = data(21)*data(51) + data(22)*data(52) + data(29)*data(48) + data(30)*data(49);
    coeffs(7, 1) = data(21)*data(39) + data(22)*data(40);
    coeffs(8, 1) = data(21)*data(48) + data(22)*data(49);
    coeffs(9, 1) = data(29)*data(36) + data(30)*data(37) + data(13)*data(57) + data(14)*data(58);
    coeffs(10, 1) = data(15) + data(13)*data(54) + data(14)*data(55);
    coeffs(11, 1) = data(13)*data(51) + data(14)*data(52) + data(29)*data(45) + data(30)*data(46);
    coeffs(12, 1) = data(13)*data(39) + data(14)*data(40) + data(21)*data(36) + data(22)*data(37);
    coeffs(13, 1) = data(13)*data(48) + data(14)*data(49) + data(21)*data(45) + data(22)*data(46);
    coeffs(14, 1) = data(13)*data(36) + data(14)*data(37);
    coeffs(15, 1) = data(13)*data(45) + data(14)*data(46);
    coeffs(16, 1) = data(5)*data(57) + data(29)*data(33) + data(6)*data(58) + data(30)*data(34);
    coeffs(17, 1) = data(7) + data(5)*data(54) + data(6)*data(55);
    coeffs(18, 1) = data(5)*data(51) + data(6)*data(52) + data(29)*data(42) + data(30)*data(43);
    coeffs(19, 1) = data(5)*data(39) + data(6)*data(40) + data(21)*data(33) + data(22)*data(34);
    coeffs(20, 1) = data(5)*data(48) + data(6)*data(49) + data(21)*data(42) + data(22)*data(43);
    coeffs(21, 1) = data(5)*data(36) + data(6)*data(37) + data(13)*data(33) + data(14)*data(34);
    coeffs(22, 1) = data(5)*data(45) + data(6)*data(46) + data(13)*data(42) + data(14)*data(43);
    coeffs(23, 1) = data(5)*data(33) + data(6)*data(34);
    coeffs(24, 1) = data(5)*data(42) + data(6)*data(43);
    coeffs(25, 1) = data(25)*data(57) + data(26)*data(58);
    coeffs(26, 1) = data(27) + data(25)*data(54) + data(26)*data(55);
    coeffs(27, 1) = data(25)*data(51) + data(26)*data(52);
    coeffs(28, 1) = data(25)*data(39) + data(26)*data(40) + data(17)*data(57) + data(18)*data(58);
    coeffs(29, 1) = data(19) + data(17)*data(54) + data(18)*data(55);
    coeffs(30, 1) = data(17)*data(51) + data(18)*data(52) + data(25)*data(48) + data(26)*data(49);
    coeffs(31, 1) = data(17)*data(39) + data(18)*data(40);
    coeffs(32, 1) = data(17)*data(48) + data(18)*data(49);
    coeffs(33, 1) = data(25)*data(36) + data(26)*data(37) + data(9)*data(57) + data(10)*data(58);
    coeffs(34, 1) = data(11) + data(9)*data(54) + data(10)*data(55);
    coeffs(35, 1) = data(9)*data(51) + data(10)*data(52) + data(25)*data(45) + data(26)*data(46);
    coeffs(36, 1) = data(9)*data(39) + data(10)*data(40) + data(17)*data(36) + data(18)*data(37);
    coeffs(37, 1) = data(9)*data(48) + data(10)*data(49) + data(17)*data(45) + data(18)*data(46);
    coeffs(38, 1) = data(9)*data(36) + data(10)*data(37);
    coeffs(39, 1) = data(9)*data(45) + data(10)*data(46);
    coeffs(40, 1) = data(1)*data(57) + data(25)*data(33) + data(2)*data(58) + data(26)*data(34);
    coeffs(41, 1) = data(3) + data(1)*data(54) + data(2)*data(55);
    coeffs(42, 1) = data(1)*data(51) + data(2)*data(52) + data(25)*data(42) + data(26)*data(43);
    coeffs(43, 1) = data(1)*data(39) + data(2)*data(40) + data(17)*data(33) + data(18)*data(34);
    coeffs(44, 1) = data(1)*data(48) + data(2)*data(49) + data(17)*data(42) + data(18)*data(43);
    coeffs(45, 1) = data(1)*data(36) + data(2)*data(37) + data(9)*data(33) + data(10)*data(34);
    coeffs(46, 1) = data(1)*data(45) + data(2)*data(46) + data(9)*data(42) + data(10)*data(43);
    coeffs(47, 1) = data(1)*data(33) + data(2)*data(34);
    coeffs(48, 1) = data(1)*data(42) + data(2)*data(43);
    coeffs(49, 1) = data(25)*data(29) + data(26)*data(30) + data(27)*data(31);
    coeffs(50, 1) = data(17)*data(29) + data(21)*data(25) + data(18)*data(30) + data(22)*data(26) + data(19)*data(31) + data(23)*data(27);
    coeffs(51, 1) = data(17)*data(21) + data(18)*data(22) + data(19)*data(23);
    coeffs(52, 1) = data(9)*data(29) + data(13)*data(25) + data(10)*data(30) + data(14)*data(26) + data(11)*data(31) + data(15)*data(27);
    coeffs(53, 1) = data(9)*data(21) + data(13)*data(17) + data(10)*data(22) + data(14)*data(18) + data(11)*data(23) + data(15)*data(19);
    coeffs(54, 1) = data(9)*data(13) + data(10)*data(14) + data(11)*data(15);
    coeffs(55, 1) = data(1)*data(29) + data(5)*data(25) + data(2)*data(30) + data(6)*data(26) + data(3)*data(31) + data(7)*data(27);
    coeffs(56, 1) = data(1)*data(21) + data(5)*data(17) + data(2)*data(22) + data(6)*data(18) + data(3)*data(23) + data(7)*data(19);
    coeffs(57, 1) = data(1)*data(13) + data(5)*data(9) + data(2)*data(14) + data(6)*data(10) + data(3)*data(15) + data(7)*data(11);
    coeffs(58, 1) = data(1)*data(5) + data(2)*data(6) + data(3)*data(7);
    coeffs(59, 1) = t2 - t3 - t4 + t14 - t15 + data(25) * data(25);
    coeffs(60, 1) = t16 + data(17)*data(25)*2.0 + data(19)*data(27)*2.0 - data(21)*data(29)*2.0 - data(22)*data(30)*2.0 - data(23)*data(31)*2.0;
    coeffs(61, 1) = t5 - t6 - t7 + t17 - t18 + data(17) * data(17);
    coeffs(62, 1) = t19 + data(9)*data(25)*2.0 + data(11)*data(27)*2.0 - data(13)*data(29)*2.0 - data(14)*data(30)*2.0 - data(15)*data(31)*2.0;
    coeffs(63, 1) = t20 + data(9)*data(17)*2.0 + data(11)*data(19)*2.0 - data(13)*data(21)*2.0 - data(14)*data(22)*2.0 - data(15)*data(23)*2.0;
    coeffs(64, 1) = t8 - t9 - t10 + t21 - t22 + data(9) * data(9);
    coeffs(65, 1) = t23 + data(1)*data(25)*2.0 + data(3)*data(27)*2.0 - data(5)*data(29)*2.0 - data(6)*data(30)*2.0 - data(7)*data(31)*2.0;
    coeffs(66, 1) = t24 + data(1)*data(17)*2.0 + data(3)*data(19)*2.0 - data(5)*data(21)*2.0 - data(6)*data(22)*2.0 - data(7)*data(23)*2.0;
    coeffs(67, 1) = t25 + data(1)*data(9)*2.0 + data(3)*data(11)*2.0 - data(5)*data(13)*2.0 - data(6)*data(14)*2.0 - data(7)*data(15)*2.0;
    coeffs(68, 1) = t11 - t12 - t13 + t26 - t27 + data(1) * data(1);
    coeffs(69, 1) = t2*data(58) - t3*data(58) - t4*data(58);
    coeffs(70, 1) = t2*data(55) - t3*data(55) - t4*data(55) - data(26)*data(27) - data(30)*data(31);
    coeffs(71, 1) = t2*data(52) - t3*data(52) - t4*data(52);
    coeffs(72, 1) = t2*data(40) - t3*data(40) - t4*data(40) + data(19)*data(27)*data(58)*2.0 - data(21)*data(29)*data(58)*2.0 - data(22)*data(30)*data(58)*2.0;
    coeffs(73, 1) = -data(18)*data(27) - data(19)*data(26) - data(22)*data(31) - data(23)*data(30) + data(19)*data(27)*data(55)*2.0 - data(21)*data(29)*data(55)*2.0 - data(22)*data(30)*data(55)*2.0;
    coeffs(74, 1) = t2*data(49) - t3*data(49) - t4*data(49) + data(19)*data(27)*data(52)*2.0 - data(21)*data(29)*data(52)*2.0 - data(22)*data(30)*data(52)*2.0;
    coeffs(75, 1) = t5*data(58) - t6*data(58) - t7*data(58) + data(19)*data(27)*data(40)*2.0 - data(21)*data(29)*data(40)*2.0 - data(22)*data(30)*data(40)*2.0;
    coeffs(76, 1) = t5*data(55) - t6*data(55) - t7*data(55) - data(18)*data(19) - data(22)*data(23);
    coeffs(77, 1) = t5*data(52) - t6*data(52) - t7*data(52) + data(19)*data(27)*data(49)*2.0 - data(21)*data(29)*data(49)*2.0 - data(22)*data(30)*data(49)*2.0;
    coeffs(78, 1) = t5*data(40) - t6*data(40) - t7*data(40);
    coeffs(79, 1) = t5*data(49) - t6*data(49) - t7*data(49);
    coeffs(80, 1) = t2*data(37) - t3*data(37) - t4*data(37) + data(11)*data(27)*data(58)*2.0 - data(13)*data(29)*data(58)*2.0 - data(14)*data(30)*data(58)*2.0;
    coeffs(81, 1) = -data(10)*data(27) - data(11)*data(26) - data(14)*data(31) - data(15)*data(30) + data(11)*data(27)*data(55)*2.0 - data(13)*data(29)*data(55)*2.0 - data(14)*data(30)*data(55)*2.0;
    coeffs(82, 1) = t2*data(46) - t3*data(46) - t4*data(46) + data(11)*data(27)*data(52)*2.0 - data(13)*data(29)*data(52)*2.0 - data(14)*data(30)*data(52)*2.0;
    coeffs(83, 1) = data(11)*data(27)*data(40)*2.0 - data(13)*data(29)*data(40)*2.0 + data(19)*data(27)*data(37)*2.0 - data(14)*data(30)*data(40)*2.0 - data(21)*data(29)*data(37)*2.0 + data(11)*data(19)*data(58)*2.0 - data(22)*data(30)*data(37)*2.0 - data(13)*data(21)*data(58)*2.0 - data(14)*data(22)*data(58)*2.0;
    coeffs(84, 1) = -data(10)*data(19) - data(11)*data(18) - data(14)*data(23) - data(15)*data(22) + data(11)*data(19)*data(55)*2.0 - data(13)*data(21)*data(55)*2.0 - data(14)*data(22)*data(55)*2.0;
    coeffs(85, 1) = data(11)*data(19)*data(52)*2.0 - data(13)*data(21)*data(52)*2.0 + data(11)*data(27)*data(49)*2.0 - data(14)*data(22)*data(52)*2.0 - data(13)*data(29)*data(49)*2.0 + data(19)*data(27)*data(46)*2.0 - data(14)*data(30)*data(49)*2.0 - data(21)*data(29)*data(46)*2.0 - data(22)*data(30)*data(46)*2.0;
    coeffs(86, 1) = t5*data(37) - t6*data(37) - t7*data(37) + data(11)*data(19)*data(40)*2.0 - data(13)*data(21)*data(40)*2.0 - data(14)*data(22)*data(40)*2.0;
    coeffs(87, 1) = t5*data(46) - t6*data(46) - t7*data(46) + data(11)*data(19)*data(49)*2.0 - data(13)*data(21)*data(49)*2.0 - data(14)*data(22)*data(49)*2.0;
    coeffs(88, 1) = t8*data(58) - t9*data(58) - t10*data(58) + data(11)*data(27)*data(37)*2.0 - data(13)*data(29)*data(37)*2.0 - data(14)*data(30)*data(37)*2.0;
    coeffs(89, 1) = t8*data(55) - t9*data(55) - t10*data(55) - data(10)*data(11) - data(14)*data(15);
    coeffs(90, 1) = t8*data(52) - t9*data(52) - t10*data(52) + data(11)*data(27)*data(46)*2.0 - data(13)*data(29)*data(46)*2.0 - data(14)*data(30)*data(46)*2.0;
    coeffs(91, 1) = t8*data(40) - t9*data(40) - t10*data(40) + data(11)*data(19)*data(37)*2.0 - data(13)*data(21)*data(37)*2.0 - data(14)*data(22)*data(37)*2.0;
    coeffs(92, 1) = t8*data(49) - t9*data(49) - t10*data(49) + data(11)*data(19)*data(46)*2.0 - data(13)*data(21)*data(46)*2.0 - data(14)*data(22)*data(46)*2.0;
    coeffs(93, 1) = t8*data(37) - t9*data(37) - t10*data(37);
    coeffs(94, 1) = t8*data(46) - t9*data(46) - t10*data(46);
    coeffs(95, 1) = t2*data(34) - t3*data(34) - t4*data(34) + data(3)*data(27)*data(58)*2.0 - data(5)*data(29)*data(58)*2.0 - data(6)*data(30)*data(58)*2.0;
    coeffs(96, 1) = -data(2)*data(27) - data(3)*data(26) - data(6)*data(31) - data(7)*data(30) + data(3)*data(27)*data(55)*2.0 - data(5)*data(29)*data(55)*2.0 - data(6)*data(30)*data(55)*2.0;
    coeffs(97, 1) = t2*data(43) - t3*data(43) - t4*data(43) + data(3)*data(27)*data(52)*2.0 - data(5)*data(29)*data(52)*2.0 - data(6)*data(30)*data(52)*2.0;
    coeffs(98, 1) = data(3)*data(27)*data(40)*2.0 - data(5)*data(29)*data(40)*2.0 - data(6)*data(30)*data(40)*2.0 + data(3)*data(19)*data(58)*2.0 + data(19)*data(27)*data(34)*2.0 - data(5)*data(21)*data(58)*2.0 - data(21)*data(29)*data(34)*2.0 - data(6)*data(22)*data(58)*2.0 - data(22)*data(30)*data(34)*2.0;
    coeffs(99, 1) = -data(2)*data(19) - data(3)*data(18) - data(6)*data(23) - data(7)*data(22) + data(3)*data(19)*data(55)*2.0 - data(5)*data(21)*data(55)*2.0 - data(6)*data(22)*data(55)*2.0;
    coeffs(100, 1) = data(3)*data(19)*data(52)*2.0 - data(5)*data(21)*data(52)*2.0 + data(3)*data(27)*data(49)*2.0 - data(6)*data(22)*data(52)*2.0 - data(5)*data(29)*data(49)*2.0 - data(6)*data(30)*data(49)*2.0 + data(19)*data(27)*data(43)*2.0 - data(21)*data(29)*data(43)*2.0 - data(22)*data(30)*data(43)*2.0;
    coeffs(101, 1) = t5*data(34) - t6*data(34) - t7*data(34) + data(3)*data(19)*data(40)*2.0 - data(5)*data(21)*data(40)*2.0 - data(6)*data(22)*data(40)*2.0;
    coeffs(102, 1) = t5*data(43) - t6*data(43) - t7*data(43) + data(3)*data(19)*data(49)*2.0 - data(5)*data(21)*data(49)*2.0 - data(6)*data(22)*data(49)*2.0;
    coeffs(103, 1) = data(3)*data(27)*data(37)*2.0 - data(5)*data(29)*data(37)*2.0 + data(3)*data(11)*data(58)*2.0 + data(11)*data(27)*data(34)*2.0 - data(6)*data(30)*data(37)*2.0 - data(5)*data(13)*data(58)*2.0 - data(13)*data(29)*data(34)*2.0 - data(6)*data(14)*data(58)*2.0 - data(14)*data(30)*data(34)*2.0;
    coeffs(104, 1) = -data(2)*data(11) - data(3)*data(10) - data(6)*data(15) - data(7)*data(14) + data(3)*data(11)*data(55)*2.0 - data(5)*data(13)*data(55)*2.0 - data(6)*data(14)*data(55)*2.0;
    coeffs(105, 1) = data(3)*data(11)*data(52)*2.0 - data(5)*data(13)*data(52)*2.0 - data(6)*data(14)*data(52)*2.0 + data(3)*data(27)*data(46)*2.0 - data(5)*data(29)*data(46)*2.0 + data(11)*data(27)*data(43)*2.0 - data(6)*data(30)*data(46)*2.0 - data(13)*data(29)*data(43)*2.0 - data(14)*data(30)*data(43)*2.0;
    coeffs(106, 1) = data(3)*data(11)*data(40)*2.0 - data(5)*data(13)*data(40)*2.0 + data(3)*data(19)*data(37)*2.0 - data(6)*data(14)*data(40)*2.0 - data(5)*data(21)*data(37)*2.0 + data(11)*data(19)*data(34)*2.0 - data(6)*data(22)*data(37)*2.0 - data(13)*data(21)*data(34)*2.0 - data(14)*data(22)*data(34)*2.0;
    coeffs(107, 1) = data(3)*data(11)*data(49)*2.0 - data(5)*data(13)*data(49)*2.0 + data(3)*data(19)*data(46)*2.0 - data(6)*data(14)*data(49)*2.0 - data(5)*data(21)*data(46)*2.0 + data(11)*data(19)*data(43)*2.0 - data(6)*data(22)*data(46)*2.0 - data(13)*data(21)*data(43)*2.0 - data(14)*data(22)*data(43)*2.0;
    coeffs(108, 1) = t8*data(34) - t9*data(34) - t10*data(34) + data(3)*data(11)*data(37)*2.0 - data(5)*data(13)*data(37)*2.0 - data(6)*data(14)*data(37)*2.0;
    coeffs(109, 1) = t8*data(43) - t9*data(43) - t10*data(43) + data(3)*data(11)*data(46)*2.0 - data(5)*data(13)*data(46)*2.0 - data(6)*data(14)*data(46)*2.0;
    coeffs(110, 1) = t11*data(58) - t12*data(58) - t13*data(58) + data(3)*data(27)*data(34)*2.0 - data(5)*data(29)*data(34)*2.0 - data(6)*data(30)*data(34)*2.0;
    coeffs(111, 1) = t11*data(55) - t12*data(55) - t13*data(55) - data(2)*data(3) - data(6)*data(7);
    coeffs(112, 1) = t11*data(52) - t12*data(52) - t13*data(52) + data(3)*data(27)*data(43)*2.0 - data(5)*data(29)*data(43)*2.0 - data(6)*data(30)*data(43)*2.0;
    coeffs(113, 1) = t11*data(40) - t12*data(40) - t13*data(40) + data(3)*data(19)*data(34)*2.0 - data(5)*data(21)*data(34)*2.0 - data(6)*data(22)*data(34)*2.0;
    coeffs(114, 1) = t11*data(49) - t12*data(49) - t13*data(49) + data(3)*data(19)*data(43)*2.0 - data(5)*data(21)*data(43)*2.0 - data(6)*data(22)*data(43)*2.0;
    coeffs(115, 1) = t11*data(37) - t12*data(37) - t13*data(37) + data(3)*data(11)*data(34)*2.0 - data(5)*data(13)*data(34)*2.0 - data(6)*data(14)*data(34)*2.0;
    coeffs(116, 1) = t11*data(46) - t12*data(46) - t13*data(46) + data(3)*data(11)*data(43)*2.0 - data(5)*data(13)*data(43)*2.0 - data(6)*data(14)*data(43)*2.0;
    coeffs(117, 1) = t11*data(34) - t12*data(34) - t13*data(34);
    coeffs(118, 1) = t11*data(43) - t12*data(43) - t13*data(43);
    coeffs(119, 1) = data(26)*data(27)*data(58) + data(30)*data(31)*data(58);
    coeffs(120, 1) = t3 - t14 + t15 + data(26)*data(27)*data(55) + data(30)*data(31)*data(55);
    coeffs(121, 1) = data(26)*data(27)*data(52) + data(30)*data(31)*data(52);
    coeffs(122, 1) = data(26)*data(27)*data(40) + data(30)*data(31)*data(40) + data(18)*data(27)*data(58) + data(19)*data(26)*data(58) + data(22)*data(31)*data(58) + data(23)*data(30)*data(58);
    coeffs(123, 1) = -t16 + data(21)*data(29)*2.0 + data(23)*data(31)*2.0 + data(18)*data(27)*data(55) + data(19)*data(26)*data(55) + data(22)*data(31)*data(55) + data(23)*data(30)*data(55);
    coeffs(124, 1) = data(18)*data(27)*data(52) + data(19)*data(26)*data(52) + data(26)*data(27)*data(49) + data(22)*data(31)*data(52) + data(23)*data(30)*data(52) + data(30)*data(31)*data(49);
    coeffs(125, 1) = data(18)*data(27)*data(40) + data(19)*data(26)*data(40) + data(22)*data(31)*data(40) + data(23)*data(30)*data(40) + data(18)*data(19)*data(58) + data(22)*data(23)*data(58);
    coeffs(126, 1) = t6 - t17 + t18 + data(18)*data(19)*data(55) + data(22)*data(23)*data(55);
    coeffs(127, 1) = data(18)*data(19)*data(52) + data(18)*data(27)*data(49) + data(19)*data(26)*data(49) + data(22)*data(23)*data(52) + data(22)*data(31)*data(49) + data(23)*data(30)*data(49);
    coeffs(128, 1) = data(18)*data(19)*data(40) + data(22)*data(23)*data(40);
    coeffs(129, 1) = data(18)*data(19)*data(49) + data(22)*data(23)*data(49);
    coeffs(130, 1) = data(26)*data(27)*data(37) + data(10)*data(27)*data(58) + data(11)*data(26)*data(58) + data(30)*data(31)*data(37) + data(14)*data(31)*data(58) + data(15)*data(30)*data(58);
    coeffs(131, 1) = -t19 + data(13)*data(29)*2.0 + data(15)*data(31)*2.0 + data(10)*data(27)*data(55) + data(11)*data(26)*data(55) + data(14)*data(31)*data(55) + data(15)*data(30)*data(55);
    coeffs(132, 1) = data(10)*data(27)*data(52) + data(11)*data(26)*data(52) + data(14)*data(31)*data(52) + data(15)*data(30)*data(52) + data(26)*data(27)*data(46) + data(30)*data(31)*data(46);
    coeffs(133, 1) = data(10)*data(27)*data(40) + data(11)*data(26)*data(40) + data(18)*data(27)*data(37) + data(19)*data(26)*data(37) + data(14)*data(31)*data(40) + data(15)*data(30)*data(40) + data(10)*data(19)*data(58) + data(11)*data(18)*data(58) + data(22)*data(31)*data(37) + data(23)*data(30)*data(37) + data(14)*data(23)*data(58) + data(15)*data(22)*data(58);
    coeffs(134, 1) = -t20 + data(13)*data(21)*2.0 + data(15)*data(23)*2.0 + data(10)*data(19)*data(55) + data(11)*data(18)*data(55) + data(14)*data(23)*data(55) + data(15)*data(22)*data(55);
    coeffs(135, 1) = data(10)*data(19)*data(52) + data(11)*data(18)*data(52) + data(10)*data(27)*data(49) + data(11)*data(26)*data(49) + data(14)*data(23)*data(52) + data(15)*data(22)*data(52) + data(18)*data(27)*data(46) + data(19)*data(26)*data(46) + data(14)*data(31)*data(49) + data(15)*data(30)*data(49) + data(22)*data(31)*data(46) + data(23)*data(30)*data(46);
    coeffs(136, 1) = data(10)*data(19)*data(40) + data(11)*data(18)*data(40) + data(18)*data(19)*data(37) + data(14)*data(23)*data(40) + data(15)*data(22)*data(40) + data(22)*data(23)*data(37);
    coeffs(137, 1) = data(10)*data(19)*data(49) + data(11)*data(18)*data(49) + data(18)*data(19)*data(46) + data(14)*data(23)*data(49) + data(15)*data(22)*data(49) + data(22)*data(23)*data(46);
    coeffs(138, 1) = data(10)*data(27)*data(37) + data(11)*data(26)*data(37) + data(10)*data(11)*data(58) + data(14)*data(31)*data(37) + data(15)*data(30)*data(37) + data(14)*data(15)*data(58);
    coeffs(139, 1) = t9 - t21 + t22 + data(10)*data(11)*data(55) + data(14)*data(15)*data(55);
    coeffs(140, 1) = data(10)*data(11)*data(52) + data(14)*data(15)*data(52) + data(10)*data(27)*data(46) + data(11)*data(26)*data(46) + data(14)*data(31)*data(46) + data(15)*data(30)*data(46);
    coeffs(141, 1) = data(10)*data(11)*data(40) + data(10)*data(19)*data(37) + data(11)*data(18)*data(37) + data(14)*data(15)*data(40) + data(14)*data(23)*data(37) + data(15)*data(22)*data(37);
    coeffs(142, 1) = data(10)*data(11)*data(49) + data(10)*data(19)*data(46) + data(11)*data(18)*data(46) + data(14)*data(15)*data(49) + data(14)*data(23)*data(46) + data(15)*data(22)*data(46);
    coeffs(143, 1) = data(10)*data(11)*data(37) + data(14)*data(15)*data(37);
    coeffs(144, 1) = data(10)*data(11)*data(46) + data(14)*data(15)*data(46);
    coeffs(145, 1) = data(2)*data(27)*data(58) + data(3)*data(26)*data(58) + data(26)*data(27)*data(34) + data(6)*data(31)*data(58) + data(7)*data(30)*data(58) + data(30)*data(31)*data(34);
    coeffs(146, 1) = -t23 + data(5)*data(29)*2.0 + data(7)*data(31)*2.0 + data(2)*data(27)*data(55) + data(3)*data(26)*data(55) + data(6)*data(31)*data(55) + data(7)*data(30)*data(55);
    coeffs(147, 1) = data(2)*data(27)*data(52) + data(3)*data(26)*data(52) + data(6)*data(31)*data(52) + data(7)*data(30)*data(52) + data(26)*data(27)*data(43) + data(30)*data(31)*data(43);
    coeffs(148, 1) = data(2)*data(27)*data(40) + data(3)*data(26)*data(40) + data(6)*data(31)*data(40) + data(7)*data(30)*data(40) + data(2)*data(19)*data(58) + data(3)*data(18)*data(58) + data(18)*data(27)*data(34) + data(19)*data(26)*data(34) + data(6)*data(23)*data(58) + data(7)*data(22)*data(58) + data(22)*data(31)*data(34) + data(23)*data(30)*data(34);
    coeffs(149, 1) = -t24 + data(5)*data(21)*2.0 + data(7)*data(23)*2.0 + data(2)*data(19)*data(55) + data(3)*data(18)*data(55) + data(6)*data(23)*data(55) + data(7)*data(22)*data(55);
    coeffs(150, 1) = data(2)*data(19)*data(52) + data(3)*data(18)*data(52) + data(2)*data(27)*data(49) + data(3)*data(26)*data(49) + data(6)*data(23)*data(52) + data(7)*data(22)*data(52) + data(6)*data(31)*data(49) + data(7)*data(30)*data(49) + data(18)*data(27)*data(43) + data(19)*data(26)*data(43) + data(22)*data(31)*data(43) + data(23)*data(30)*data(43);
    coeffs(151, 1) = data(2)*data(19)*data(40) + data(3)*data(18)*data(40) + data(6)*data(23)*data(40) + data(7)*data(22)*data(40) + data(18)*data(19)*data(34) + data(22)*data(23)*data(34);
    coeffs(152, 1) = data(2)*data(19)*data(49) + data(3)*data(18)*data(49) + data(6)*data(23)*data(49) + data(7)*data(22)*data(49) + data(18)*data(19)*data(43) + data(22)*data(23)*data(43);
    coeffs(153, 1) = data(2)*data(27)*data(37) + data(3)*data(26)*data(37) + data(2)*data(11)*data(58) + data(3)*data(10)*data(58) + data(10)*data(27)*data(34) + data(11)*data(26)*data(34) + data(6)*data(31)*data(37) + data(7)*data(30)*data(37) + data(6)*data(15)*data(58) + data(7)*data(14)*data(58) + data(14)*data(31)*data(34) + data(15)*data(30)*data(34);
    coeffs(154, 1) = -t25 + data(5)*data(13)*2.0 + data(7)*data(15)*2.0 + data(2)*data(11)*data(55) + data(3)*data(10)*data(55) + data(6)*data(15)*data(55) + data(7)*data(14)*data(55);
    coeffs(155, 1) = data(2)*data(11)*data(52) + data(3)*data(10)*data(52) + data(6)*data(15)*data(52) + data(7)*data(14)*data(52) + data(2)*data(27)*data(46) + data(3)*data(26)*data(46) + data(10)*data(27)*data(43) + data(11)*data(26)*data(43) + data(6)*data(31)*data(46) + data(7)*data(30)*data(46) + data(14)*data(31)*data(43) + data(15)*data(30)*data(43);
    coeffs(156, 1) = data(2)*data(11)*data(40) + data(3)*data(10)*data(40) + data(2)*data(19)*data(37) + data(3)*data(18)*data(37) + data(6)*data(15)*data(40) + data(7)*data(14)*data(40) + data(10)*data(19)*data(34) + data(11)*data(18)*data(34) + data(6)*data(23)*data(37) + data(7)*data(22)*data(37) + data(14)*data(23)*data(34) + data(15)*data(22)*data(34);
    coeffs(157, 1) = data(2)*data(11)*data(49) + data(3)*data(10)*data(49) + data(2)*data(19)*data(46) + data(3)*data(18)*data(46) + data(6)*data(15)*data(49) + data(7)*data(14)*data(49) + data(10)*data(19)*data(43) + data(11)*data(18)*data(43) + data(6)*data(23)*data(46) + data(7)*data(22)*data(46) + data(14)*data(23)*data(43) + data(15)*data(22)*data(43);
    coeffs(158, 1) = data(2)*data(11)*data(37) + data(3)*data(10)*data(37) + data(10)*data(11)*data(34) + data(6)*data(15)*data(37) + data(7)*data(14)*data(37) + data(14)*data(15)*data(34);
    coeffs(159, 1) = data(2)*data(11)*data(46) + data(3)*data(10)*data(46) + data(10)*data(11)*data(43) + data(6)*data(15)*data(46) + data(7)*data(14)*data(46) + data(14)*data(15)*data(43);
    coeffs(160, 1) = data(2)*data(3)*data(58) + data(2)*data(27)*data(34) + data(3)*data(26)*data(34) + data(6)*data(7)*data(58) + data(6)*data(31)*data(34) + data(7)*data(30)*data(34);
    coeffs(161, 1) = t12 - t26 + t27 + data(2)*data(3)*data(55) + data(6)*data(7)*data(55);
    coeffs(162, 1) = data(2)*data(3)*data(52) + data(6)*data(7)*data(52) + data(2)*data(27)*data(43) + data(3)*data(26)*data(43) + data(6)*data(31)*data(43) + data(7)*data(30)*data(43);
    coeffs(163, 1) = data(2)*data(3)*data(40) + data(6)*data(7)*data(40) + data(2)*data(19)*data(34) + data(3)*data(18)*data(34) + data(6)*data(23)*data(34) + data(7)*data(22)*data(34);
    coeffs(164, 1) = data(2)*data(3)*data(49) + data(6)*data(7)*data(49) + data(2)*data(19)*data(43) + data(3)*data(18)*data(43) + data(6)*data(23)*data(43) + data(7)*data(22)*data(43);
    coeffs(165, 1) = data(2)*data(3)*data(37) + data(2)*data(11)*data(34) + data(3)*data(10)*data(34) + data(6)*data(7)*data(37) + data(6)*data(15)*data(34) + data(7)*data(14)*data(34);
    coeffs(166, 1) = data(2)*data(3)*data(46) + data(2)*data(11)*data(43) + data(3)*data(10)*data(43) + data(6)*data(7)*data(46) + data(6)*data(15)*data(43) + data(7)*data(14)*data(43);
    coeffs(167, 1) = data(2)*data(3)*data(34) + data(6)*data(7)*data(34);
    coeffs(168, 1) = data(2)*data(3)*data(43) + data(6)*data(7)*data(43);
    coeffs(169, 1) = data(25)*data(27)*data(58) + data(29)*data(31)*data(58);
    coeffs(170, 1) = -data(25)*data(26) - data(29)*data(30) + data(25)*data(27)*data(55) + data(29)*data(31)*data(55);
    coeffs(171, 1) = data(25)*data(27)*data(52) + data(29)*data(31)*data(52);
    coeffs(172, 1) = data(25)*data(27)*data(40) + data(29)*data(31)*data(40) + data(17)*data(27)*data(58) + data(19)*data(25)*data(58) + data(21)*data(31)*data(58) + data(23)*data(29)*data(58);
    coeffs(173, 1) = -data(17)*data(26) - data(18)*data(25) - data(21)*data(30) - data(22)*data(29) + data(17)*data(27)*data(55) + data(19)*data(25)*data(55) + data(21)*data(31)*data(55) + data(23)*data(29)*data(55);
    coeffs(174, 1) = data(17)*data(27)*data(52) + data(19)*data(25)*data(52) + data(25)*data(27)*data(49) + data(21)*data(31)*data(52) + data(23)*data(29)*data(52) + data(29)*data(31)*data(49);
    coeffs(175, 1) = data(17)*data(27)*data(40) + data(19)*data(25)*data(40) + data(21)*data(31)*data(40) + data(23)*data(29)*data(40) + data(17)*data(19)*data(58) + data(21)*data(23)*data(58);
    coeffs(176, 1) = -data(17)*data(18) - data(21)*data(22) + data(17)*data(19)*data(55) + data(21)*data(23)*data(55);
    coeffs(177, 1) = data(17)*data(19)*data(52) + data(17)*data(27)*data(49) + data(19)*data(25)*data(49) + data(21)*data(23)*data(52) + data(21)*data(31)*data(49) + data(23)*data(29)*data(49);
    coeffs(178, 1) = data(17)*data(19)*data(40) + data(21)*data(23)*data(40);
    coeffs(179, 1) = data(17)*data(19)*data(49) + data(21)*data(23)*data(49);
    coeffs(180, 1) = data(25)*data(27)*data(37) + data(9)*data(27)*data(58) + data(11)*data(25)*data(58) + data(29)*data(31)*data(37) + data(13)*data(31)*data(58) + data(15)*data(29)*data(58);
    coeffs(181, 1) = -data(9)*data(26) - data(10)*data(25) - data(13)*data(30) - data(14)*data(29) + data(9)*data(27)*data(55) + data(11)*data(25)*data(55) + data(13)*data(31)*data(55) + data(15)*data(29)*data(55);
    coeffs(182, 1) = data(9)*data(27)*data(52) + data(11)*data(25)*data(52) + data(13)*data(31)*data(52) + data(15)*data(29)*data(52) + data(25)*data(27)*data(46) + data(29)*data(31)*data(46);
    coeffs(183, 1) = data(9)*data(27)*data(40) + data(11)*data(25)*data(40) + data(17)*data(27)*data(37) + data(19)*data(25)*data(37) + data(13)*data(31)*data(40) + data(15)*data(29)*data(40) + data(9)*data(19)*data(58) + data(11)*data(17)*data(58) + data(21)*data(31)*data(37) + data(23)*data(29)*data(37) + data(13)*data(23)*data(58) + data(15)*data(21)*data(58);
    coeffs(184, 1) = -data(9)*data(18) - data(10)*data(17) - data(13)*data(22) - data(14)*data(21) + data(9)*data(19)*data(55) + data(11)*data(17)*data(55) + data(13)*data(23)*data(55) + data(15)*data(21)*data(55);
    coeffs(185, 1) = data(9)*data(19)*data(52) + data(11)*data(17)*data(52) + data(9)*data(27)*data(49) + data(11)*data(25)*data(49) + data(13)*data(23)*data(52) + data(15)*data(21)*data(52) + data(17)*data(27)*data(46) + data(19)*data(25)*data(46) + data(13)*data(31)*data(49) + data(15)*data(29)*data(49) + data(21)*data(31)*data(46) + data(23)*data(29)*data(46);
    coeffs(186, 1) = data(9)*data(19)*data(40) + data(11)*data(17)*data(40) + data(17)*data(19)*data(37) + data(13)*data(23)*data(40) + data(15)*data(21)*data(40) + data(21)*data(23)*data(37);
    coeffs(187, 1) = data(9)*data(19)*data(49) + data(11)*data(17)*data(49) + data(17)*data(19)*data(46) + data(13)*data(23)*data(49) + data(15)*data(21)*data(49) + data(21)*data(23)*data(46);
    coeffs(188, 1) = data(9)*data(27)*data(37) + data(11)*data(25)*data(37) + data(9)*data(11)*data(58) + data(13)*data(31)*data(37) + data(15)*data(29)*data(37) + data(13)*data(15)*data(58);
    coeffs(189, 1) = -data(9)*data(10) - data(13)*data(14) + data(9)*data(11)*data(55) + data(13)*data(15)*data(55);
    coeffs(190, 1) = data(9)*data(11)*data(52) + data(13)*data(15)*data(52) + data(9)*data(27)*data(46) + data(11)*data(25)*data(46) + data(13)*data(31)*data(46) + data(15)*data(29)*data(46);
    coeffs(191, 1) = data(9)*data(11)*data(40) + data(9)*data(19)*data(37) + data(11)*data(17)*data(37) + data(13)*data(15)*data(40) + data(13)*data(23)*data(37) + data(15)*data(21)*data(37);
    coeffs(192, 1) = data(9)*data(11)*data(49) + data(9)*data(19)*data(46) + data(11)*data(17)*data(46) + data(13)*data(15)*data(49) + data(13)*data(23)*data(46) + data(15)*data(21)*data(46);
    coeffs(193, 1) = data(9)*data(11)*data(37) + data(13)*data(15)*data(37);
    coeffs(194, 1) = data(9)*data(11)*data(46) + data(13)*data(15)*data(46);
    coeffs(195, 1) = data(1)*data(27)*data(58) + data(3)*data(25)*data(58) + data(25)*data(27)*data(34) + data(5)*data(31)*data(58) + data(7)*data(29)*data(58) + data(29)*data(31)*data(34);
    coeffs(196, 1) = -data(1)*data(26) - data(2)*data(25) - data(5)*data(30) - data(6)*data(29) + data(1)*data(27)*data(55) + data(3)*data(25)*data(55) + data(5)*data(31)*data(55) + data(7)*data(29)*data(55);
    coeffs(197, 1) = data(1)*data(27)*data(52) + data(3)*data(25)*data(52) + data(5)*data(31)*data(52) + data(7)*data(29)*data(52) + data(25)*data(27)*data(43) + data(29)*data(31)*data(43);
    coeffs(198, 1) = data(1)*data(27)*data(40) + data(3)*data(25)*data(40) + data(5)*data(31)*data(40) + data(7)*data(29)*data(40) + data(1)*data(19)*data(58) + data(3)*data(17)*data(58) + data(17)*data(27)*data(34) + data(19)*data(25)*data(34) + data(5)*data(23)*data(58) + data(7)*data(21)*data(58) + data(21)*data(31)*data(34) + data(23)*data(29)*data(34);
    coeffs(199, 1) = -data(1)*data(18) - data(2)*data(17) - data(5)*data(22) - data(6)*data(21) + data(1)*data(19)*data(55) + data(3)*data(17)*data(55) + data(5)*data(23)*data(55) + data(7)*data(21)*data(55);
    coeffs(200, 1) = data(1)*data(19)*data(52) + data(3)*data(17)*data(52) + data(1)*data(27)*data(49) + data(3)*data(25)*data(49) + data(5)*data(23)*data(52) + data(7)*data(21)*data(52) + data(5)*data(31)*data(49) + data(7)*data(29)*data(49) + data(17)*data(27)*data(43) + data(19)*data(25)*data(43) + data(21)*data(31)*data(43) + data(23)*data(29)*data(43);
    coeffs(201, 1) = data(1)*data(19)*data(40) + data(3)*data(17)*data(40) + data(5)*data(23)*data(40) + data(7)*data(21)*data(40) + data(17)*data(19)*data(34) + data(21)*data(23)*data(34);
    coeffs(202, 1) = data(1)*data(19)*data(49) + data(3)*data(17)*data(49) + data(5)*data(23)*data(49) + data(7)*data(21)*data(49) + data(17)*data(19)*data(43) + data(21)*data(23)*data(43);
    coeffs(203, 1) = data(1)*data(27)*data(37) + data(3)*data(25)*data(37) + data(1)*data(11)*data(58) + data(3)*data(9)*data(58) + data(9)*data(27)*data(34) + data(11)*data(25)*data(34) + data(5)*data(31)*data(37) + data(7)*data(29)*data(37) + data(5)*data(15)*data(58) + data(7)*data(13)*data(58) + data(13)*data(31)*data(34) + data(15)*data(29)*data(34);
    coeffs(204, 1) = -data(1)*data(10) - data(2)*data(9) - data(5)*data(14) - data(6)*data(13) + data(1)*data(11)*data(55) + data(3)*data(9)*data(55) + data(5)*data(15)*data(55) + data(7)*data(13)*data(55);
    coeffs(205, 1) = data(1)*data(11)*data(52) + data(3)*data(9)*data(52) + data(5)*data(15)*data(52) + data(7)*data(13)*data(52) + data(1)*data(27)*data(46) + data(3)*data(25)*data(46) + data(9)*data(27)*data(43) + data(11)*data(25)*data(43) + data(5)*data(31)*data(46) + data(7)*data(29)*data(46) + data(13)*data(31)*data(43) + data(15)*data(29)*data(43);
    coeffs(206, 1) = data(1)*data(11)*data(40) + data(3)*data(9)*data(40) + data(1)*data(19)*data(37) + data(3)*data(17)*data(37) + data(5)*data(15)*data(40) + data(7)*data(13)*data(40) + data(9)*data(19)*data(34) + data(11)*data(17)*data(34) + data(5)*data(23)*data(37) + data(7)*data(21)*data(37) + data(13)*data(23)*data(34) + data(15)*data(21)*data(34);
    coeffs(207, 1) = data(1)*data(11)*data(49) + data(3)*data(9)*data(49) + data(1)*data(19)*data(46) + data(3)*data(17)*data(46) + data(5)*data(15)*data(49) + data(7)*data(13)*data(49) + data(9)*data(19)*data(43) + data(11)*data(17)*data(43) + data(5)*data(23)*data(46) + data(7)*data(21)*data(46) + data(13)*data(23)*data(43) + data(15)*data(21)*data(43);
    coeffs(208, 1) = data(1)*data(11)*data(37) + data(3)*data(9)*data(37) + data(9)*data(11)*data(34) + data(5)*data(15)*data(37) + data(7)*data(13)*data(37) + data(13)*data(15)*data(34);
    coeffs(209, 1) = data(1)*data(11)*data(46) + data(3)*data(9)*data(46) + data(9)*data(11)*data(43) + data(5)*data(15)*data(46) + data(7)*data(13)*data(46) + data(13)*data(15)*data(43);
    coeffs(210, 1) = data(1)*data(3)*data(58) + data(1)*data(27)*data(34) + data(3)*data(25)*data(34) + data(5)*data(7)*data(58) + data(5)*data(31)*data(34) + data(7)*data(29)*data(34);
    coeffs(211, 1) = -data(1)*data(2) - data(5)*data(6) + data(1)*data(3)*data(55) + data(5)*data(7)*data(55);
    coeffs(212, 1) = data(1)*data(3)*data(52) + data(5)*data(7)*data(52) + data(1)*data(27)*data(43) + data(3)*data(25)*data(43) + data(5)*data(31)*data(43) + data(7)*data(29)*data(43);
    coeffs(213, 1) = data(1)*data(3)*data(40) + data(5)*data(7)*data(40) + data(1)*data(19)*data(34) + data(3)*data(17)*data(34) + data(5)*data(23)*data(34) + data(7)*data(21)*data(34);
    coeffs(214, 1) = data(1)*data(3)*data(49) + data(5)*data(7)*data(49) + data(1)*data(19)*data(43) + data(3)*data(17)*data(43) + data(5)*data(23)*data(43) + data(7)*data(21)*data(43);
    coeffs(215, 1) = data(1)*data(3)*data(37) + data(1)*data(11)*data(34) + data(3)*data(9)*data(34) + data(5)*data(7)*data(37) + data(5)*data(15)*data(34) + data(7)*data(13)*data(34);
    coeffs(216, 1) = data(1)*data(3)*data(46) + data(1)*data(11)*data(43) + data(3)*data(9)*data(43) + data(5)*data(7)*data(46) + data(5)*data(15)*data(43) + data(7)*data(13)*data(43);
    coeffs(217, 1) = data(1)*data(3)*data(34) + data(5)*data(7)*data(34);
    coeffs(218, 1) = data(1)*data(3)*data(43) + data(5)*data(7)*data(43);
    coeffs(219, 1) = t2*data(57) - t4*data(57) + data(29)*data(30)*data(58);
    coeffs(220, 1) = t2*data(54) - t4*data(54) - data(25)*data(27) + data(29)*data(30)*data(55);
    coeffs(221, 1) = t2*data(51) - t4*data(51) + data(29)*data(30)*data(52);
    coeffs(222, 1) = t2*data(39) - t4*data(39) + data(29)*data(30)*data(40) + data(19)*data(27)*data(57)*2.0 + data(21)*data(30)*data(58) + data(22)*data(29)*data(58) - data(22)*data(30)*data(57)*2.0;
    coeffs(223, 1) = -data(17)*data(27) - data(19)*data(25) + data(19)*data(27)*data(54)*2.0 + data(21)*data(30)*data(55) + data(22)*data(29)*data(55) - data(22)*data(30)*data(54)*2.0;
    coeffs(224, 1) = t2*data(48) - t4*data(48) + data(19)*data(27)*data(51)*2.0 + data(21)*data(30)*data(52) + data(22)*data(29)*data(52) - data(22)*data(30)*data(51)*2.0 + data(29)*data(30)*data(49);
    coeffs(225, 1) = t5*data(57) - t7*data(57) + data(19)*data(27)*data(39)*2.0 + data(21)*data(30)*data(40) + data(22)*data(29)*data(40) - data(22)*data(30)*data(39)*2.0 + data(21)*data(22)*data(58);
    coeffs(226, 1) = t5*data(54) - t7*data(54) - data(17)*data(19) + data(21)*data(22)*data(55);
    coeffs(227, 1) = t5*data(51) - t7*data(51) + data(19)*data(27)*data(48)*2.0 + data(21)*data(22)*data(52) + data(21)*data(30)*data(49) + data(22)*data(29)*data(49) - data(22)*data(30)*data(48)*2.0;
    coeffs(228, 1) = t5*data(39) - t7*data(39) + data(21)*data(22)*data(40);
    coeffs(229, 1) = t5*data(48) - t7*data(48) + data(21)*data(22)*data(49);
    coeffs(230, 1) = t2*data(36) - t4*data(36) + data(11)*data(27)*data(57)*2.0 + data(29)*data(30)*data(37) + data(13)*data(30)*data(58) + data(14)*data(29)*data(58) - data(14)*data(30)*data(57)*2.0;
    coeffs(231, 1) = -data(9)*data(27) - data(11)*data(25) + data(11)*data(27)*data(54)*2.0 + data(13)*data(30)*data(55) + data(14)*data(29)*data(55) - data(14)*data(30)*data(54)*2.0;
    coeffs(232, 1) = t2*data(45) - t4*data(45) + data(11)*data(27)*data(51)*2.0 + data(13)*data(30)*data(52) + data(14)*data(29)*data(52) - data(14)*data(30)*data(51)*2.0 + data(29)*data(30)*data(46);
    coeffs(233, 1) = data(11)*data(27)*data(39)*2.0 + data(19)*data(27)*data(36)*2.0 + data(13)*data(30)*data(40) + data(14)*data(29)*data(40) - data(14)*data(30)*data(39)*2.0 + data(11)*data(19)*data(57)*2.0 + data(21)*data(30)*data(37) + data(22)*data(29)*data(37) - data(22)*data(30)*data(36)*2.0 + data(13)*data(22)*data(58) + data(14)*data(21)*data(58) - data(14)*data(22)*data(57)*2.0;
    coeffs(234, 1) = -data(9)*data(19) - data(11)*data(17) + data(11)*data(19)*data(54)*2.0 + data(13)*data(22)*data(55) + data(14)*data(21)*data(55) - data(14)*data(22)*data(54)*2.0;
    coeffs(235, 1) = data(11)*data(19)*data(51)*2.0 + data(11)*data(27)*data(48)*2.0 + data(13)*data(22)*data(52) + data(14)*data(21)*data(52) - data(14)*data(22)*data(51)*2.0 + data(19)*data(27)*data(45)*2.0 + data(13)*data(30)*data(49) + data(14)*data(29)*data(49) - data(14)*data(30)*data(48)*2.0 + data(21)*data(30)*data(46) + data(22)*data(29)*data(46) - data(22)*data(30)*data(45)*2.0;
    coeffs(236, 1) = t5*data(36) - t7*data(36) + data(11)*data(19)*data(39)*2.0 + data(13)*data(22)*data(40) + data(14)*data(21)*data(40) - data(14)*data(22)*data(39)*2.0 + data(21)*data(22)*data(37);
    coeffs(237, 1) = t5*data(45) - t7*data(45) + data(11)*data(19)*data(48)*2.0 + data(13)*data(22)*data(49) + data(14)*data(21)*data(49) - data(14)*data(22)*data(48)*2.0 + data(21)*data(22)*data(46);
    coeffs(238, 1) = t8*data(57) - t10*data(57) + data(11)*data(27)*data(36)*2.0 + data(13)*data(30)*data(37) + data(14)*data(29)*data(37) - data(14)*data(30)*data(36)*2.0 + data(13)*data(14)*data(58);
    coeffs(239, 1) = t8*data(54) - t10*data(54) - data(9)*data(11) + data(13)*data(14)*data(55);
    coeffs(240, 1) = t8*data(51) - t10*data(51) + data(13)*data(14)*data(52) + data(11)*data(27)*data(45)*2.0 + data(13)*data(30)*data(46) + data(14)*data(29)*data(46) - data(14)*data(30)*data(45)*2.0;
    coeffs(241, 1) = t8*data(39) - t10*data(39) + data(11)*data(19)*data(36)*2.0 + data(13)*data(14)*data(40) + data(13)*data(22)*data(37) + data(14)*data(21)*data(37) - data(14)*data(22)*data(36)*2.0;
    coeffs(242, 1) = t8*data(48) - t10*data(48) + data(11)*data(19)*data(45)*2.0 + data(13)*data(14)*data(49) + data(13)*data(22)*data(46) + data(14)*data(21)*data(46) - data(14)*data(22)*data(45)*2.0;
    coeffs(243, 1) = t8*data(36) - t10*data(36) + data(13)*data(14)*data(37);
    coeffs(244, 1) = t8*data(45) - t10*data(45) + data(13)*data(14)*data(46);
    coeffs(245, 1) = t2*data(33) - t4*data(33) + data(3)*data(27)*data(57)*2.0 + data(5)*data(30)*data(58) + data(6)*data(29)*data(58) - data(6)*data(30)*data(57)*2.0 + data(29)*data(30)*data(34);
    coeffs(246, 1) = -data(1)*data(27) - data(3)*data(25) + data(3)*data(27)*data(54)*2.0 + data(5)*data(30)*data(55) + data(6)*data(29)*data(55) - data(6)*data(30)*data(54)*2.0;
    coeffs(247, 1) = t2*data(42) - t4*data(42) + data(3)*data(27)*data(51)*2.0 + data(5)*data(30)*data(52) + data(6)*data(29)*data(52) - data(6)*data(30)*data(51)*2.0 + data(29)*data(30)*data(43);
    coeffs(248, 1) = data(3)*data(27)*data(39)*2.0 + data(5)*data(30)*data(40) + data(6)*data(29)*data(40) - data(6)*data(30)*data(39)*2.0 + data(3)*data(19)*data(57)*2.0 + data(19)*data(27)*data(33)*2.0 + data(5)*data(22)*data(58) + data(6)*data(21)*data(58) - data(6)*data(22)*data(57)*2.0 + data(21)*data(30)*data(34) + data(22)*data(29)*data(34) - data(22)*data(30)*data(33)*2.0;
    coeffs(249, 1) = -data(1)*data(19) - data(3)*data(17) + data(3)*data(19)*data(54)*2.0 + data(5)*data(22)*data(55) + data(6)*data(21)*data(55) - data(6)*data(22)*data(54)*2.0;
    coeffs(250, 1) = data(3)*data(19)*data(51)*2.0 + data(3)*data(27)*data(48)*2.0 + data(5)*data(22)*data(52) + data(6)*data(21)*data(52) - data(6)*data(22)*data(51)*2.0 + data(5)*data(30)*data(49) + data(6)*data(29)*data(49) - data(6)*data(30)*data(48)*2.0 + data(19)*data(27)*data(42)*2.0 + data(21)*data(30)*data(43) + data(22)*data(29)*data(43) - data(22)*data(30)*data(42)*2.0;
    coeffs(251, 1) = t5*data(33) - t7*data(33) + data(3)*data(19)*data(39)*2.0 + data(5)*data(22)*data(40) + data(6)*data(21)*data(40) - data(6)*data(22)*data(39)*2.0 + data(21)*data(22)*data(34);
    coeffs(252, 1) = t5*data(42) - t7*data(42) + data(3)*data(19)*data(48)*2.0 + data(5)*data(22)*data(49) + data(6)*data(21)*data(49) - data(6)*data(22)*data(48)*2.0 + data(21)*data(22)*data(43);
    coeffs(253, 1) = data(3)*data(27)*data(36)*2.0 + data(3)*data(11)*data(57)*2.0 + data(11)*data(27)*data(33)*2.0 + data(5)*data(30)*data(37) + data(6)*data(29)*data(37) - data(6)*data(30)*data(36)*2.0 + data(5)*data(14)*data(58) + data(6)*data(13)*data(58) - data(6)*data(14)*data(57)*2.0 + data(13)*data(30)*data(34) + data(14)*data(29)*data(34) - data(14)*data(30)*data(33)*2.0;
    coeffs(254, 1) = -data(1)*data(11) - data(3)*data(9) + data(3)*data(11)*data(54)*2.0 + data(5)*data(14)*data(55) + data(6)*data(13)*data(55) - data(6)*data(14)*data(54)*2.0;
    coeffs(255, 1) = data(3)*data(11)*data(51)*2.0 + data(5)*data(14)*data(52) + data(6)*data(13)*data(52) - data(6)*data(14)*data(51)*2.0 + data(3)*data(27)*data(45)*2.0 + data(11)*data(27)*data(42)*2.0 + data(5)*data(30)*data(46) + data(6)*data(29)*data(46) - data(6)*data(30)*data(45)*2.0 + data(13)*data(30)*data(43) + data(14)*data(29)*data(43) - data(14)*data(30)*data(42)*2.0;
    coeffs(256, 1) = data(3)*data(11)*data(39)*2.0 + data(3)*data(19)*data(36)*2.0 + data(5)*data(14)*data(40) + data(6)*data(13)*data(40) - data(6)*data(14)*data(39)*2.0 + data(11)*data(19)*data(33)*2.0 + data(5)*data(22)*data(37) + data(6)*data(21)*data(37) - data(6)*data(22)*data(36)*2.0 + data(13)*data(22)*data(34) + data(14)*data(21)*data(34) - data(14)*data(22)*data(33)*2.0;
    coeffs(257, 1) = data(3)*data(11)*data(48)*2.0 + data(3)*data(19)*data(45)*2.0 + data(5)*data(14)*data(49) + data(6)*data(13)*data(49) - data(6)*data(14)*data(48)*2.0 + data(11)*data(19)*data(42)*2.0 + data(5)*data(22)*data(46) + data(6)*data(21)*data(46) - data(6)*data(22)*data(45)*2.0 + data(13)*data(22)*data(43) + data(14)*data(21)*data(43) - data(14)*data(22)*data(42)*2.0;
    coeffs(258, 1) = t8*data(33) - t10*data(33) + data(3)*data(11)*data(36)*2.0 + data(5)*data(14)*data(37) + data(6)*data(13)*data(37) - data(6)*data(14)*data(36)*2.0 + data(13)*data(14)*data(34);
    coeffs(259, 1) = t8*data(42) - t10*data(42) + data(3)*data(11)*data(45)*2.0 + data(5)*data(14)*data(46) + data(6)*data(13)*data(46) - data(6)*data(14)*data(45)*2.0 + data(13)*data(14)*data(43);
    coeffs(260, 1) = t11*data(57) - t13*data(57) + data(3)*data(27)*data(33)*2.0 + data(5)*data(6)*data(58) + data(5)*data(30)*data(34) + data(6)*data(29)*data(34) - data(6)*data(30)*data(33)*2.0;
    coeffs(261, 1) = t11*data(54) - t13*data(54) - data(1)*data(3) + data(5)*data(6)*data(55);
    coeffs(262, 1) = t11*data(51) - t13*data(51) + data(5)*data(6)*data(52) + data(3)*data(27)*data(42)*2.0 + data(5)*data(30)*data(43) + data(6)*data(29)*data(43) - data(6)*data(30)*data(42)*2.0;
    coeffs(263, 1) = t11*data(39) - t13*data(39) + data(5)*data(6)*data(40) + data(3)*data(19)*data(33)*2.0 + data(5)*data(22)*data(34) + data(6)*data(21)*data(34) - data(6)*data(22)*data(33)*2.0;
    coeffs(264, 1) = t11*data(48) - t13*data(48) + data(5)*data(6)*data(49) + data(3)*data(19)*data(42)*2.0 + data(5)*data(22)*data(43) + data(6)*data(21)*data(43) - data(6)*data(22)*data(42)*2.0;
    coeffs(265, 1) = t11*data(36) - t13*data(36) + data(3)*data(11)*data(33)*2.0 + data(5)*data(6)*data(37) + data(5)*data(14)*data(34) + data(6)*data(13)*data(34) - data(6)*data(14)*data(33)*2.0;
    coeffs(266, 1) = t11*data(45) - t13*data(45) + data(3)*data(11)*data(42)*2.0 + data(5)*data(6)*data(46) + data(5)*data(14)*data(43) + data(6)*data(13)*data(43) - data(6)*data(14)*data(42)*2.0;
    coeffs(267, 1) = t11*data(33) - t13*data(33) + data(5)*data(6)*data(34);
    coeffs(268, 1) = t11*data(42) - t13*data(42) + data(5)*data(6)*data(43);
    coeffs(269, 1) = data(26)*data(27)*data(57) + data(30)*data(31)*data(57);
    coeffs(270, 1) = -data(25)*data(26) - data(29)*data(30) + data(26)*data(27)*data(54) + data(30)*data(31)*data(54);
    coeffs(271, 1) = data(26)*data(27)*data(51) + data(30)*data(31)*data(51);
    coeffs(272, 1) = data(26)*data(27)*data(39) + data(30)*data(31)*data(39) + data(18)*data(27)*data(57) + data(19)*data(26)*data(57) + data(22)*data(31)*data(57) + data(23)*data(30)*data(57);
    coeffs(273, 1) = -data(17)*data(26) - data(18)*data(25) - data(21)*data(30) - data(22)*data(29) + data(18)*data(27)*data(54) + data(19)*data(26)*data(54) + data(22)*data(31)*data(54) + data(23)*data(30)*data(54);
    coeffs(274, 1) = data(18)*data(27)*data(51) + data(19)*data(26)*data(51) + data(26)*data(27)*data(48) + data(22)*data(31)*data(51) + data(23)*data(30)*data(51) + data(30)*data(31)*data(48);
    coeffs(275, 1) = data(18)*data(27)*data(39) + data(19)*data(26)*data(39) + data(22)*data(31)*data(39) + data(23)*data(30)*data(39) + data(18)*data(19)*data(57) + data(22)*data(23)*data(57);
    coeffs(276, 1) = -data(17)*data(18) - data(21)*data(22) + data(18)*data(19)*data(54) + data(22)*data(23)*data(54);
    coeffs(277, 1) = data(18)*data(19)*data(51) + data(18)*data(27)*data(48) + data(19)*data(26)*data(48) + data(22)*data(23)*data(51) + data(22)*data(31)*data(48) + data(23)*data(30)*data(48);
    coeffs(278, 1) = data(18)*data(19)*data(39) + data(22)*data(23)*data(39);
    coeffs(279, 1) = data(18)*data(19)*data(48) + data(22)*data(23)*data(48);
    coeffs(280, 1) = data(26)*data(27)*data(36) + data(10)*data(27)*data(57) + data(11)*data(26)*data(57) + data(30)*data(31)*data(36) + data(14)*data(31)*data(57) + data(15)*data(30)*data(57);
    coeffs(281, 1) = -data(9)*data(26) - data(10)*data(25) - data(13)*data(30) - data(14)*data(29) + data(10)*data(27)*data(54) + data(11)*data(26)*data(54) + data(14)*data(31)*data(54) + data(15)*data(30)*data(54);
    coeffs(282, 1) = data(10)*data(27)*data(51) + data(11)*data(26)*data(51) + data(14)*data(31)*data(51) + data(15)*data(30)*data(51) + data(26)*data(27)*data(45) + data(30)*data(31)*data(45);
    coeffs(283, 1) = data(10)*data(27)*data(39) + data(11)*data(26)*data(39) + data(18)*data(27)*data(36) + data(19)*data(26)*data(36) + data(14)*data(31)*data(39) + data(15)*data(30)*data(39) + data(10)*data(19)*data(57) + data(11)*data(18)*data(57) + data(22)*data(31)*data(36) + data(23)*data(30)*data(36) + data(14)*data(23)*data(57) + data(15)*data(22)*data(57);
    coeffs(284, 1) = -data(9)*data(18) - data(10)*data(17) - data(13)*data(22) - data(14)*data(21) + data(10)*data(19)*data(54) + data(11)*data(18)*data(54) + data(14)*data(23)*data(54) + data(15)*data(22)*data(54);
    coeffs(285, 1) = data(10)*data(19)*data(51) + data(11)*data(18)*data(51) + data(10)*data(27)*data(48) + data(11)*data(26)*data(48) + data(14)*data(23)*data(51) + data(15)*data(22)*data(51) + data(18)*data(27)*data(45) + data(19)*data(26)*data(45) + data(14)*data(31)*data(48) + data(15)*data(30)*data(48) + data(22)*data(31)*data(45) + data(23)*data(30)*data(45);
    coeffs(286, 1) = data(10)*data(19)*data(39) + data(11)*data(18)*data(39) + data(18)*data(19)*data(36) + data(14)*data(23)*data(39) + data(15)*data(22)*data(39) + data(22)*data(23)*data(36);
    coeffs(287, 1) = data(10)*data(19)*data(48) + data(11)*data(18)*data(48) + data(18)*data(19)*data(45) + data(14)*data(23)*data(48) + data(15)*data(22)*data(48) + data(22)*data(23)*data(45);
    coeffs(288, 1) = data(10)*data(27)*data(36) + data(11)*data(26)*data(36) + data(10)*data(11)*data(57) + data(14)*data(31)*data(36) + data(15)*data(30)*data(36) + data(14)*data(15)*data(57);
    coeffs(289, 1) = -data(9)*data(10) - data(13)*data(14) + data(10)*data(11)*data(54) + data(14)*data(15)*data(54);
    coeffs(290, 1) = data(10)*data(11)*data(51) + data(14)*data(15)*data(51) + data(10)*data(27)*data(45) + data(11)*data(26)*data(45) + data(14)*data(31)*data(45) + data(15)*data(30)*data(45);
    coeffs(291, 1) = data(10)*data(11)*data(39) + data(10)*data(19)*data(36) + data(11)*data(18)*data(36) + data(14)*data(15)*data(39) + data(14)*data(23)*data(36) + data(15)*data(22)*data(36);
    coeffs(292, 1) = data(10)*data(11)*data(48) + data(10)*data(19)*data(45) + data(11)*data(18)*data(45) + data(14)*data(15)*data(48) + data(14)*data(23)*data(45) + data(15)*data(22)*data(45);
    coeffs(293, 1) = data(10)*data(11)*data(36) + data(14)*data(15)*data(36);
    coeffs(294, 1) = data(10)*data(11)*data(45) + data(14)*data(15)*data(45);
    coeffs(295, 1) = data(2)*data(27)*data(57) + data(3)*data(26)*data(57) + data(26)*data(27)*data(33) + data(6)*data(31)*data(57) + data(7)*data(30)*data(57) + data(30)*data(31)*data(33);
    coeffs(296, 1) = -data(1)*data(26) - data(2)*data(25) - data(5)*data(30) - data(6)*data(29) + data(2)*data(27)*data(54) + data(3)*data(26)*data(54) + data(6)*data(31)*data(54) + data(7)*data(30)*data(54);
    coeffs(297, 1) = data(2)*data(27)*data(51) + data(3)*data(26)*data(51) + data(6)*data(31)*data(51) + data(7)*data(30)*data(51) + data(26)*data(27)*data(42) + data(30)*data(31)*data(42);
    coeffs(298, 1) = data(2)*data(27)*data(39) + data(3)*data(26)*data(39) + data(6)*data(31)*data(39) + data(7)*data(30)*data(39) + data(2)*data(19)*data(57) + data(3)*data(18)*data(57) + data(18)*data(27)*data(33) + data(19)*data(26)*data(33) + data(6)*data(23)*data(57) + data(7)*data(22)*data(57) + data(22)*data(31)*data(33) + data(23)*data(30)*data(33);
    coeffs(299, 1) = -data(1)*data(18) - data(2)*data(17) - data(5)*data(22) - data(6)*data(21) + data(2)*data(19)*data(54) + data(3)*data(18)*data(54) + data(6)*data(23)*data(54) + data(7)*data(22)*data(54);
    coeffs(300, 1) = data(2)*data(19)*data(51) + data(3)*data(18)*data(51) + data(2)*data(27)*data(48) + data(3)*data(26)*data(48) + data(6)*data(23)*data(51) + data(7)*data(22)*data(51) + data(6)*data(31)*data(48) + data(7)*data(30)*data(48) + data(18)*data(27)*data(42) + data(19)*data(26)*data(42) + data(22)*data(31)*data(42) + data(23)*data(30)*data(42);
    coeffs(301, 1) = data(2)*data(19)*data(39) + data(3)*data(18)*data(39) + data(6)*data(23)*data(39) + data(7)*data(22)*data(39) + data(18)*data(19)*data(33) + data(22)*data(23)*data(33);
    coeffs(302, 1) = data(2)*data(19)*data(48) + data(3)*data(18)*data(48) + data(6)*data(23)*data(48) + data(7)*data(22)*data(48) + data(18)*data(19)*data(42) + data(22)*data(23)*data(42);
    coeffs(303, 1) = data(2)*data(27)*data(36) + data(3)*data(26)*data(36) + data(2)*data(11)*data(57) + data(3)*data(10)*data(57) + data(10)*data(27)*data(33) + data(11)*data(26)*data(33) + data(6)*data(31)*data(36) + data(7)*data(30)*data(36) + data(6)*data(15)*data(57) + data(7)*data(14)*data(57) + data(14)*data(31)*data(33) + data(15)*data(30)*data(33);
    coeffs(304, 1) = -data(1)*data(10) - data(2)*data(9) - data(5)*data(14) - data(6)*data(13) + data(2)*data(11)*data(54) + data(3)*data(10)*data(54) + data(6)*data(15)*data(54) + data(7)*data(14)*data(54);
    coeffs(305, 1) = data(2)*data(11)*data(51) + data(3)*data(10)*data(51) + data(6)*data(15)*data(51) + data(7)*data(14)*data(51) + data(2)*data(27)*data(45) + data(3)*data(26)*data(45) + data(10)*data(27)*data(42) + data(11)*data(26)*data(42) + data(6)*data(31)*data(45) + data(7)*data(30)*data(45) + data(14)*data(31)*data(42) + data(15)*data(30)*data(42);
    coeffs(306, 1) = data(2)*data(11)*data(39) + data(3)*data(10)*data(39) + data(2)*data(19)*data(36) + data(3)*data(18)*data(36) + data(6)*data(15)*data(39) + data(7)*data(14)*data(39) + data(10)*data(19)*data(33) + data(11)*data(18)*data(33) + data(6)*data(23)*data(36) + data(7)*data(22)*data(36) + data(14)*data(23)*data(33) + data(15)*data(22)*data(33);
    coeffs(307, 1) = data(2)*data(11)*data(48) + data(3)*data(10)*data(48) + data(2)*data(19)*data(45) + data(3)*data(18)*data(45) + data(6)*data(15)*data(48) + data(7)*data(14)*data(48) + data(10)*data(19)*data(42) + data(11)*data(18)*data(42) + data(6)*data(23)*data(45) + data(7)*data(22)*data(45) + data(14)*data(23)*data(42) + data(15)*data(22)*data(42);
    coeffs(308, 1) = data(2)*data(11)*data(36) + data(3)*data(10)*data(36) + data(10)*data(11)*data(33) + data(6)*data(15)*data(36) + data(7)*data(14)*data(36) + data(14)*data(15)*data(33);
    coeffs(309, 1) = data(2)*data(11)*data(45) + data(3)*data(10)*data(45) + data(10)*data(11)*data(42) + data(6)*data(15)*data(45) + data(7)*data(14)*data(45) + data(14)*data(15)*data(42);
    coeffs(310, 1) = data(2)*data(3)*data(57) + data(2)*data(27)*data(33) + data(3)*data(26)*data(33) + data(6)*data(7)*data(57) + data(6)*data(31)*data(33) + data(7)*data(30)*data(33);
    coeffs(311, 1) = -data(1)*data(2) - data(5)*data(6) + data(2)*data(3)*data(54) + data(6)*data(7)*data(54);
    coeffs(312, 1) = data(2)*data(3)*data(51) + data(6)*data(7)*data(51) + data(2)*data(27)*data(42) + data(3)*data(26)*data(42) + data(6)*data(31)*data(42) + data(7)*data(30)*data(42);
    coeffs(313, 1) = data(2)*data(3)*data(39) + data(6)*data(7)*data(39) + data(2)*data(19)*data(33) + data(3)*data(18)*data(33) + data(6)*data(23)*data(33) + data(7)*data(22)*data(33);
    coeffs(314, 1) = data(2)*data(3)*data(48) + data(6)*data(7)*data(48) + data(2)*data(19)*data(42) + data(3)*data(18)*data(42) + data(6)*data(23)*data(42) + data(7)*data(22)*data(42);
    coeffs(315, 1) = data(2)*data(3)*data(36) + data(2)*data(11)*data(33) + data(3)*data(10)*data(33) + data(6)*data(7)*data(36) + data(6)*data(15)*data(33) + data(7)*data(14)*data(33);
    coeffs(316, 1) = data(2)*data(3)*data(45) + data(2)*data(11)*data(42) + data(3)*data(10)*data(42) + data(6)*data(7)*data(45) + data(6)*data(15)*data(42) + data(7)*data(14)*data(42);
    coeffs(317, 1) = data(2)*data(3)*data(33) + data(6)*data(7)*data(33);
    coeffs(318, 1) = data(2)*data(3)*data(42) + data(6)*data(7)*data(42);
    coeffs(319, 1) = -data(59) - data(57)*data(61) - data(58)*data(62) + 1.0;
    coeffs(320, 1) = -data(56) - data(63) - data(54)*data(61) - data(55)*data(62);
    coeffs(321, 1) = -data(53) + data(60) - data(51)*data(61) - data(52)*data(62);
    coeffs(322, 1) = -data(41) - data(39)*data(61) - data(40)*data(62);
    coeffs(323, 1) = -data(50) - data(48)*data(61) - data(49)*data(62);
    coeffs(324, 1) = -data(38) - data(36)*data(61) - data(37)*data(62);
    coeffs(325, 1) = -data(47) - data(45)*data(61) - data(46)*data(62);
    coeffs(326, 1) = -data(35) - data(33)*data(61) - data(34)*data(62);
    coeffs(327, 1) = -data(44) - data(42)*data(61) - data(43)*data(62);
    // way too lazy to change all the number ;-)
    coeffs1 = coeffs.bottomRows<327>().col(1);

    std::vector<int> C_ind_r = { 3,7,16,25,26,27,28,29,30,2,3,6,7,14,16,23,25,26,27,28,29,30,39,2,3,6,7,14,16,23,25,26,27,28,29,30,36,39,2,6,14,23,26,27,28,29,30,36,1,3,5,7,12,16,21,25,26,27,28,29,30,38,1,2,5,6,12,14,21,23,26,27,28,29,30,35,36,1,3,5,7,12,16,21,25,26,27,28,29,30,33,38,1,2,5,6,12,14,21,23,26,27,28,29,30,33,35,1,5,12,21,26,27,28,29,30,33,3,7,15,24,26,27,28,29,30,2,3,6,7,13,15,22,24,26,27,28,29,30,39,2,3,6,7,13,15,22,24,26,27,28,29,30,36,39,2,6,13,22,26,27,28,29,30,36,1,3,5,7,11,15,20,24,26,27,28,29,30,38,1,2,3,5,6,7,11,13,15,20,22,24,26,27,28,29,30,35,38,39,0,3,4,7,10,16,19,25,26,27,28,29,30,37,0,2,4,6,10,14,19,23,26,27,28,29,30,34,36,0,1,3,4,5,7,10,12,16,19,21,25,26,27,28,29,30,32,37,38,0,1,2,4,5,6,10,12,14,19,21,23,26,27,28,29,30,32,34,35,3,7,9,18,26,27,28,29,30,2,3,6,7,9,18,26,27,28,29,30,39,2,6,9,18,26,27,28,29,30,36,1,3,5,7,9,18,26,27,28,29,30,38,1,2,5,6,9,18,26,27,28,29,30,35,0,3,4,7,8,15,17,24,26,27,28,29,30,37,0,2,3,4,6,7,8,13,15,17,22,24,26,27,28,29,30,34,37,39,0,3,4,7,10,16,19,25,26,27,28,29,30,31,37,0,2,4,6,10,14,19,23,26,27,28,29,30,31,34,0,3,4,7,9,18,26,27,28,29,30,37,0,2,4,6,9,18,26,27,28,29,30,34,0,1,4,5,10,12,19,21,26,27,28,29,30,32,33,1,5,9,18,26,27,28,29,30,33,1,3,5,7,11,15,20,24,26,27,28,29,30,33,38,1,2,5,6,11,13,20,22,26,27,28,29,30,35,36,1,2,5,6,11,13,20,22,26,27,28,29,30,33,35,1,5,11,20,26,27,28,29,30,33,1,2,3,5,6,7,12,14,16,21,23,25,26,27,28,29,30,35,38,39,0,4,8,17,26,27,28,29,30,31,0,4,10,19,26,27,28,29,30,31,0,1,4,5,10,12,19,21,26,27,28,29,30,31,32,0,4,9,18,26,27,28,29,30,31,0,1,4,5,9,18,26,27,28,29,30,32,0,3,4,7,8,15,17,24,26,27,28,29,30,31,37,0,1,3,4,5,7,8,11,15,17,20,24,26,27,28,29,30,32,37,38,0,2,4,6,8,13,17,22,26,27,28,29,30,31,34,0,2,4,6,8,13,17,22,26,27,28,29,30,34,36,0,1,2,4,5,6,8,11,13,17,20,22,26,27,28,29,30,32,34,35,0,1,4,5,8,11,17,20,26,27,28,29,30,31,32,0,1,4,5,8,11,17,20,26,27,28,29,30,32,33,0,2,3,4,6,7,10,14,16,19,23,25,26,27,28,29,30,34,37,39
    };
    std::vector<int> C_ind_c = { 0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,28,28,28,28,28,28,28,28,28,28,28,28,29,29,29,29,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,31,31,31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,35,35,35,35,35,35,35,35,35,35,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,37,37,37,37,37,37,37,37,37,37,38,38,38,38,38,38,38,38,38,38,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,40,40,40,40,40,40,40,40,40,40,41,41,41,41,41,41,41,41,41,41,41,41,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49
    };

    std::vector<int> coeffs_ind_r = { 23,47,57,67,117,167,217,267,317,23,21,47,45,57,56,67,66,115,165,215,265,315,326,21,14,45,38,56,53,66,63,108,158,208,258,308,326,324,14,38,53,63,93,143,193,243,293,324,23,19,47,43,57,55,67,65,113,163,213,263,313,326,14,12,38,36,53,52,63,62,91,141,191,241,291,324,322,19,7,43,31,55,50,65,60,101,151,201,251,301,326,322,12,7,36,31,52,50,62,60,86,136,186,236,286,324,322,7,31,50,60,78,128,178,228,278,322,22,46,57,67,116,166,216,266,316,22,20,46,44,57,56,67,66,114,164,214,264,314,325,20,13,44,37,56,53,66,63,107,157,207,257,307,325,323,13,37,53,63,92,142,192,242,292,323,22,18,46,42,57,55,67,65,112,162,212,262,312,325,20,18,11,44,42,35,56,55,52,66,65,62,105,155,205,255,305,325,323,321,23,17,47,41,57,54,67,64,111,161,211,261,311,326,14,10,38,34,53,51,63,61,89,139,189,239,289,324,320,19,17,5,43,41,29,55,54,49,65,64,59,99,149,199,249,299,326,322,320,12,10,5,36,34,29,52,51,49,62,61,59,84,134,184,234,284,324,322,320,16,40,57,67,110,160,210,260,310,16,9,40,33,56,66,103,153,203,253,303,319,9,33,53,63,88,138,188,238,288,319,16,4,40,28,55,65,98,148,198,248,298,319,9,4,33,28,52,62,83,133,183,233,283,319,22,15,46,39,57,54,67,64,109,159,209,259,309,325,20,15,8,44,39,32,56,54,51,66,64,61,102,152,202,252,302,325,323,318,17,2,41,26,54,48,64,58,96,146,196,246,296,326,320,10,2,34,26,51,48,61,58,81,131,181,231,281,324,320,16,1,40,25,54,64,95,145,195,245,295,319,9,1,33,25,51,61,80,130,180,230,280,319,7,5,31,29,50,49,60,59,76,126,176,226,276,322,320,4,28,50,60,75,125,175,225,275,319,18,6,42,30,55,50,65,60,100,150,200,250,300,325,321,13,11,37,35,53,52,63,62,90,140,190,240,290,323,321,11,6,35,30,52,50,62,60,85,135,185,235,285,323,321,6,30,50,60,77,127,177,227,277,321,21,19,12,45,43,36,56,55,52,66,65,62,106,156,206,256,306,326,324,322,0,24,48,58,68,118,168,218,268,318,2,26,48,58,70,120,170,220,270,320,5,2,29,26,49,48,59,58,73,123,173,223,273,322,320,1,25,48,58,69,119,169,219,269,319,4,1,28,25,49,59,72,122,172,222,272,319,15,0,39,24,54,48,64,58,94,144,194,244,294,325,318,18,15,3,42,39,27,55,54,49,65,64,59,97,147,197,247,297,325,321,318,8,0,32,24,51,48,61,58,79,129,179,229,279,323,318,13,8,37,32,53,51,63,61,87,137,187,237,287,323,318,11,8,3,35,32,27,52,51,49,62,61,59,82,132,182,232,282,323,321,318,3,0,27,24,49,48,59,58,71,121,171,221,271,321,318,6,3,30,27,50,49,60,59,74,124,174,224,274,321,318,21,17,10,45,41,34,56,54,51,66,64,61,104,154,204,254,304,326,324,320
    };
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(40,50);
    for (int i = 0; i < C_ind_r.size(); ++i)
        C(C_ind_r[i], C_ind_c[i]) = coeffs1(coeffs_ind_r[i], 0);

    std::vector<int> b_ind_r = { 0  ,   1  ,   2 ,    3 ,    4   ,  5  ,   6 };
    std::vector<int> b_ind_c = { 30 ,   31 ,   32  ,  33,    34,    35 ,   36 };

    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(7,37);
    for (int i = 0; i < b_ind_r.size(); ++i)
        b(b_ind_r[i], b_ind_c[i]) = -1.0;

    Eigen::MatrixXd C0 = Eigen::MatrixXd::Zero(40, 37);
    C0 = C.leftCols<37>();
    Eigen::MatrixXd	C1 = Eigen::MatrixXd::Zero(40, 13);
    C1 = C.rightCols<13>();

    Eigen::MatrixXd alpha = C0.transpose().fullPivLu().solve(b.transpose());
    Eigen::MatrixXd RR = Eigen::MatrixXd::Zero(20,13);
    RR.block<7, 13>(0, 0) = alpha.transpose() * C1;
    RR.block<13, 13>(7, 0) = Eigen::MatrixXd::Identity(13, 13);
    std::vector<int> AM_ind = { 17  ,   9 ,    0  ,  11  ,   1  ,  13  ,   2  ,  16  ,   3  ,   4 ,   18  ,   5  ,   6 };
    Eigen::MatrixXd AM = Eigen::MatrixXd::Zero(13,13);
    for (int i = 0; i < AM_ind.size(); ++i)
        AM.row(i) = RR.row(AM_ind[i]);

    Eigen::EigenSolver<Eigen::MatrixXd> eig;
    eig.compute(AM);

    Eigen::ArrayXcd D = eig.eigenvalues();
    Eigen::ArrayXXcd V = eig.eigenvectors();
    V = (V / V.row(0).replicate(13, 1)).eval();
    sols.row(0) = V.row(5);
    sols.row(1) = V.row(7);
    sols.row(2) = D.transpose();
    sols.row(3) = V.row(1);
    sols.row(4) = V.row(3);
}


bool FourPointsPoseFocalLengthRadialDistortion(
    const Eigen::Matrix<double, 3, 4>& feature_vectors,
    const Eigen::Matrix4d& world_points,
    std::vector<Eigen::Matrix3d>& rotations,
    std::vector<Eigen::Vector3d>& translations,
    std::vector<double>& radial_distortions,
    std::vector<double>& focal_lengths)
{
    Eigen::Vector4d d = feature_vectors.topRows<2>().array().pow(2).colwise().sum();

    Eigen::Vector3d t0 = world_points.topRows<3>().rowwise().mean();
    Eigen::Matrix<double, 3, 4> t0_mat;
    t0_mat.col(0) = t0; t0_mat.col(1) = t0; t0_mat.col(2) = t0; t0_mat.col(3) = t0;
    Eigen::Matrix<double, 4, 4> U = world_points;
    U.topRows<3>() -= t0_mat;

    Eigen::JacobiSVD < Eigen::MatrixXd > svd(U.topRows<3>(), Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::Matrix3d R0 = svd.matrixU();

    if (sgn(R0.determinant()) < 0.0)
        R0.col(0) *= -1;
    R0.transposeInPlace();

    U.topRows<3>() = R0 * U.topRows<3>();
    double scale = U.topRows<3>().array().pow(2).colwise().sum().sqrt().mean();

    U.topRows<3>() /= scale;

    // rescale image points
    Eigen::Matrix<double, 3, 4> u = feature_vectors;
    double f0 = u.topRows<2>().array().pow(2).colwise().sum().sqrt().mean();
    u.topRows<2>() /= f0;

    double k0 = d.array().mean();
    d /= k0;

    Eigen::Matrix<double, 5, 8> M = Eigen::MatrixXd::Zero(5,8);
    M.row(0).leftCols<4>() = U.col(0);
    M.row(1).rightCols<4>() = U.col(0);
    for (int k = 1; k < 4; ++k)
    {
        M.row(k + 1).leftCols<4>()  = u.row(1).col(k)*U.col(k).transpose();
        M.row(k + 1).rightCols<4>() = -u.row(0).col(k)*U.col(k).transpose();
    }

    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(5,1);
    b.topRows<2>() = u.col(0).topRows<2>();
    Eigen::HouseholderQR<Eigen::MatrixXd> qr;
    qr.compute(M.transpose());

    Eigen::Matrix<double, 8, 8> Q = qr.householderQ();
    Eigen::Matrix<double, 8, 5> R = qr.matrixQR().triangularView<Eigen::Upper>();

    Eigen::Matrix<double, 8, 4> N;
    N.setZero();
    N.leftCols<3>() = Q.rightCols<3>();

    // this random rotation is supposed to make the solver more stable
    Eigen::Vector3d rot_vec(rng.RandDouble(-0.5,0.5),
                            rng.RandDouble(-0.5,0.5),
                            rng.RandDouble(-0.5,0.5));
    Eigen::AngleAxisd random_rot(rot_vec.norm(), rot_vec);
    N.leftCols<3>() *= random_rot.toRotationMatrix();
    Eigen::Matrix<double, 8, 1> x0 = Q.leftCols<5>()* (R.topRows<5>().transpose().inverse() * b);
    N.rightCols<1>() = x0;


    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 3);
    Eigen::Matrix<double, 3, 4> UN1 = U.rightCols<3>().transpose()*N.topRows<4>();
    Eigen::Matrix<double, 3, 4> UN2 = U.rightCols<3>().transpose()*N.bottomRows<4>();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 9);

    B.topLeftCorner<3, 3>() = UN1.leftCols<3>();
    B.bottomLeftCorner<3, 3>() = UN2.leftCols<3>();
    B.topRightCorner<3, 1>() = UN1.rightCols<1>();
    B.bottomRightCorner<3, 1>() = UN2.rightCols<1>();

    B.block<3, 4>(0, 3) = d.bottomRows<3>().transpose().replicate(4, 1).transpose().cwiseProduct(UN1);
    B.block<3, 4>(3, 3) = d.bottomRows<3>().transpose().replicate(4, 1).transpose().cwiseProduct(UN2);

    B.col(7).topRows<3>()    = -u.row(0).rightCols<3>().transpose().cwiseProduct(U.row(2).rightCols<3>().transpose());
    B.col(7).bottomRows<3>() = -u.row(1).rightCols<3>().transpose().cwiseProduct(U.row(2).rightCols<3>().transpose());

    // fill these guys
    Eigen::Matrix3d Utmp;
    Utmp.row(0) = U.row(0).rightCols<3>();
    Utmp.row(1) = U.row(1).rightCols<3>();
    Utmp.row(2) = U.row(3).rightCols<3>();

    Eigen::Matrix3d u1temp;
    u1temp.row(0) = u.row(0).rightCols<3>();
    u1temp.row(1) = u.row(0).rightCols<3>();
    u1temp.row(2) = u.row(0).rightCols<3>();
    Eigen::Matrix3d u2temp;
    u2temp.row(0) = u.row(1).rightCols<3>();
    u2temp.row(1) = u.row(1).rightCols<3>();
    u2temp.row(2) = u.row(1).rightCols<3>();

    C.block<3, 3>(0, 0) = Utmp.transpose().cwiseProduct(u1temp.transpose());
    C.block<3, 3>(3, 0) = Utmp.transpose().cwiseProduct(u2temp.transpose());

    Eigen::MatrixXd D = C.jacobiSvd(Eigen::ComputeThinV | Eigen::ComputeThinU).solve(B);

    Eigen::Map<Eigen::RowVectorXd> N_(N.data(), N.size());
    Eigen::Map<Eigen::RowVectorXd> D_(D.data(), D.size());

    Eigen::Matrix<double, 64, 1> data;
    data(0) = 0.0; // using this to keep matlab indices, just an index offset of 1
    data.block<32, 1>(1, 0) = N_;
    data.block<27, 1>(33, 0) = D_;
    data(60,0) = d(0,0);
    data.bottomRows<3>() = U.topRows<3>().col(0);

    Eigen::Matrix<std::complex<double>, 5, 13> sols;
    Solve(data, sols);
    std::vector<int> valid_sols;

    for (int k = 0; k < 13; ++k)
    {
        if (sols(0, k).imag() < -1e-6 ||
            sols(0, k).imag() > 1e-6)
            continue;
        else
            valid_sols.push_back(k);
    }
    const int nr_valid_sols = valid_sols.size();

    Eigen::MatrixXd sols0 = Eigen::MatrixXd::Zero(5, nr_valid_sols);

    for (int k = 0; k < nr_valid_sols; ++k)
        sols0.col(k) = sols.col(valid_sols[k]).real();
    rotations.resize(nr_valid_sols);
    translations.resize(nr_valid_sols);
    radial_distortions.resize(nr_valid_sols);
    focal_lengths.resize(nr_valid_sols);

    for (int i = 0; i < nr_valid_sols; ++i)
    {
        Eigen::Vector3d alpha = sols0.topRows<3>().col(i);
        double k = sols0(3, i);
        double P33 = sols0(4, i);
        Eigen::Vector4d alpha1;
        alpha1.topRows<3>() = alpha;
        alpha1(3) = 1.0;
        Eigen::VectorXd P12_ = N*alpha1;
        Eigen::Map<Eigen::MatrixXd> P12(P12_.data(), 4, 2);

        Eigen::VectorXd tmp = Eigen::VectorXd::Zero(9);
        tmp.topRows<3>() = alpha;
        tmp(3) = k*alpha(0);
        tmp(4) = k*alpha(1);
        tmp(5) = k*alpha(2);
        tmp(6) = k;
        tmp(7) = P33;
        tmp(8) = 1.0;

        Eigen::Vector3d P3_124 = D * tmp;

        Eigen::Vector4d P3;
        P3(0) = P3_124(0);
        P3(1) = P3_124(1);
        P3(2) = P33;
        P3(3) = P3_124(2);

        Eigen::Matrix<double, 3, 4> P;
        P.topRows<2>()    = P12.transpose();
        P.bottomRows<1>() = P3;
        P /= P.bottomRows<1>().leftCols<3>().norm();
        double f = P.topRows<1>().leftCols<3>().norm();
        Eigen::Matrix3d K = Eigen::Matrix3d::Identity();
        K(0, 0) = 1./f;
        K(1, 1) = 1./f;

        focal_lengths[i] = f;
        radial_distortions[i] = k;

        Eigen::MatrixXd Rt = K * P;

        if (Rt.topLeftCorner<3, 3>().determinant() < 0.0)
            Rt *= -1.0;

        focal_lengths[i] = focal_lengths[i] * f0;
        radial_distortions[i] = radial_distortions[i] / k0;
        // scale radial distortion
        // radial_distortions[i] *= (focal_lengths[i]*focal_lengths[i]);

        Rt.col(3) = Rt.col(3)*scale - Rt.topLeftCorner<3, 3>()*R0*t0;
        Rt.topLeftCorner<3, 3>() = Rt.topLeftCorner<3, 3>() * R0;

        rotations[i] = Rt.topLeftCorner<3, 3>();
        translations[i] = Rt.col(3);
    }


    return nr_valid_sols > 0;
}


}
