/**
 * @file utils.h
 * @brief Some utilities functions
 * @date Oct 07, 2014
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#ifndef AKAZE_UTILS_H_
#define AKAZE_UTILS_H_

/* ************************************************************************* */
#include "cimg/CImg.h"
#include "src/AKAZE.h"

// System
#include <vector>
#include <iostream>
#include <iomanip>

namespace libAKAZE {

/* ************************************************************************* */
// Stringify common types such as int, double and others.
template <typename T> inline std::string to_string(const T& x) {
  std::stringstream oss;
  oss << x;
  return oss.str();
}

// Stringify and format integral types as follows:
// to_formatted_string(  1, 2) produces string:  '01'
// to_formatted_string(  5, 2) produces string:  '05'
// to_formatted_string( 19, 2) produces string:  '19'
// to_formatted_string( 19, 3) produces string: '019'
template <typename Integer>
inline std::string to_formatted_string(Integer x, int num_digits) {
  std::stringstream oss;
  oss << std::setfill('0') << std::setw(num_digits) << x;
  return oss.str();
}

/* ************************************************************************* */

// Converts to a CImg type and adjusts the scale for image displaying.
void ConvertEigenToCImg(const RowMatrixXf& mat,
                        cimg_library::CImg<float>& cimg);

/* ************************************************************************* */

// Converts a CImg object to a grayscale floating point Eigen matrix.
void ConvertCImgToEigen(const cimg_library::CImg<float>& image,
                        RowMatrixXf& eigen_image);

/* ************************************************************************* */
// This function matches the descriptors from floating point AKAZE methods using
// L2 distance and the nearest neighbor distance ratio (Lowes ratio)
void match_features(const std::vector<Eigen::VectorXf>& desc1,
                    const std::vector<Eigen::VectorXf>& desc2,
                    const double ratio,
                    std::vector<std::pair<int, int> >& matches);
// This function matches the descriptors from binary AKAZE (MLDB) methods using
// Hamming distance and the nearest neighbor distance ratio (Lowes ratio)
void match_features(const std::vector<libAKAZE::BinaryVectorX>& desc1,
                    const std::vector<libAKAZE::BinaryVectorX>& desc2,
                    const double ratio,
                    std::vector<std::pair<int, int> >& matches);

/// This function computes the set of inliers given a ground truth homography
/// @param kpts1 Keypoints from first image
/// @param kpts2 Keypoints from second image
/// @param matches Vector of putative matches
/// @param inliers Vector of inliers
/// @param H Ground truth homography matrix 3x3
/// @param min_error The minimum pixelic error to accept an inlier
void compute_inliers_homography(
    const std::vector<libAKAZE::AKAZEKeypoint>& kpts1,
    const std::vector<libAKAZE::AKAZEKeypoint>& kpts2,
    const std::vector<std::pair<int, int> >& matches,
    std::vector<std::pair<int, int> >& inliers, const Eigen::Matrix3f& H,
    float min_error);

/// This function draws the list of detected keypoints
void draw_keypoints(cimg_library::CImg<float>& img,
                    const std::vector<libAKAZE::AKAZEKeypoint>& kpts);

void draw_matches(cimg_library::CImg<float>& image1,
                  cimg_library::CImg<float>& image2,
                  const std::vector<libAKAZE::AKAZEKeypoint>& kpts1,
                  const std::vector<libAKAZE::AKAZEKeypoint>& kpts2,
                  const std::vector<std::pair<int, int> >& matches,
                  cimg_library::CImg<float>& matched_image);

/// This function saves the interest points to a regular ASCII file
/// @note The format is compatible with Mikolajczyk and Schmid evaluation
/// @param outFile Name of the output file where the points will be stored
/// @param kpts Vector of points of interest
/// @param desc Matrix that contains the extracted descriptors
/// @param save_desc Set to 1 if we want to save the descriptors
int save_keypoints(const std::string& outFile,
                   const std::vector<libAKAZE::AKAZEKeypoint>& kpts,
                   const std::vector<Eigen::VectorXf>& desc,
                   bool save_desc);

int save_keypoints(const std::string& outFile,
                   const std::vector<libAKAZE::AKAZEKeypoint>& kpts,
                   const std::vector<libAKAZE::BinaryVectorX>& desc,
                   bool save_desc);

/// Function for reading the ground truth homography from a txt file
bool read_homography(const std::string& hFile, Eigen::Matrix3f& H1toN);

/// This function shows the possible command line configuration options
void show_input_options_help(int example);

}  // namespace libAKAZE
#endif  // AKAZE_UTILS_H_
