//=============================================================================
//
// akaze_match.cpp
// Authors: Pablo F. Alcantarilla (1), Jesus Nuevo (2)
// Institutions: Toshiba Research Europe Ltd (1)
//               TrueVision Solutions (2)
// Date: 07/10/2014
// Email: pablofdezalc@gmail.com
//
// AKAZE Features Copyright 2014, Pablo F. Alcantarilla, Jesus Nuevo
// All Rights Reserved
// See LICENSE for the license information
//=============================================================================

/**
 * @file akaze_match.cpp
 * @brief Main program for matching two images with AKAZE features
 * @date Oct 07, 2014
 * @author Pablo F. Alcantarilla
 */

#include <string>
#include <utility>
#include <vector>

#include "cimg/CImg.h"
#include "src/AKAZE.h"
#include "src/AKAZEConfig.h"
#include "src/utils.h"
#include "timer/timer.hpp"

/* ************************************************************************* */
// Image matching options
const float MIN_H_ERROR =
    2.50f;                   ///< Maximum error in pixels to accept an inlier
const float DRATIO = 0.80f;  ///< NNDR Matching value

/* ************************************************************************* */
/**
 * @brief This function parses the command line arguments for setting A-KAZE
 * parameters
 * and image matching between two input images
 * @param options Structure that contains A-KAZE settings
 * @param img_path1 Path for the first input image
 * @param img_path2 Path for the second input image
 * @param homography_path Path for the file that contains the ground truth
 * homography
 */
int parse_input_options(libAKAZE::AKAZEOptions& options, std::string& img_path1,
                        std::string& img_path2, std::string& homography_path,
                        int argc, char* argv[]);

/* ************************************************************************* */
int main(int argc, char* argv[]) {
  // Variables
  libAKAZE::AKAZEOptions options;
  std::string img_path1, img_path2, homography_path;
  float ratio = 0.0;
  int nkpts1 = 0, nkpts2 = 0, nmatches = 0, ninliers = 0, noutliers = 0;

  std::vector<libAKAZE::AKAZEKeypoint> kpts1, kpts2;
  std::vector<std::pair<int, int> > dmatches;
  libAKAZE::AKAZEDescriptors desc1, desc2;
  Eigen::Matrix3f HG;

  // Variables for measuring computation times
  double takaze = 0.0, tmatch = 0.0;

  // Parse the input command line options
  if (parse_input_options(options, img_path1, img_path2, homography_path, argc,
                          argv)) {
    return -1;
  }

  // Read image 1 and if necessary convert to grayscale.
  cimg_library::CImg<float> img1(img_path1.c_str());
  cimg_library::CImg<float> img2(img_path2.c_str());

  // Read ground truth homography file
  if (libAKAZE::read_homography(homography_path, HG) == false) {
    std::cout << "Invalid homography file!" << std::endl;
    return -1;
  }

  // Convert the images to float
  RowMatrixXf img1_32, img2_32;
  libAKAZE::ConvertCImgToEigen(img1, img1_32);
  img1_32 /= 255.0;
  libAKAZE::ConvertCImgToEigen(img2, img2_32);
  img2_32 /= 255.0;

  // Create the first AKAZE object
  options.img_width = img1.width();
  options.img_height = img1.height();
  libAKAZE::AKAZE evolution1(options);

  // Create the second HKAZE object
  options.img_width = img2.width();
  options.img_height = img2.height();
  libAKAZE::AKAZE evolution2(options);

  timer::Timer timer;

  // Create the nonlinear scale space
  // and perform feature detection and description for image 1
  evolution1.Create_Nonlinear_Scale_Space(img1_32);
  evolution1.Feature_Detection(kpts1);
  evolution1.Compute_Descriptors(kpts1, desc1);
  evolution2.Create_Nonlinear_Scale_Space(img2_32);
  evolution2.Feature_Detection(kpts2);
  evolution2.Compute_Descriptors(kpts2, desc2);

  takaze = timer.elapsedMs();

  nkpts1 = kpts1.size();
  nkpts2 = kpts2.size();

  // Matching Descriptors!!
  std::vector<std::pair<int, int> > matches, inliers;
  timer.reset();
  if (options.descriptor < libAKAZE::MLDB_UPRIGHT) {
    libAKAZE::match_features(desc1.float_descriptor, desc2.float_descriptor, DRATIO,
                   matches);
  }
  // Binary descriptor, use Hamming distance
  else {
    libAKAZE::match_features(desc1.binary_descriptor, desc2.binary_descriptor, DRATIO,
                   matches);
  }

  tmatch = timer.elapsedMs();

  // Compute Inliers!!
  libAKAZE::compute_inliers_homography(kpts1, kpts2, matches, inliers, HG, MIN_H_ERROR);

  // Compute the inliers statistics
  nmatches = matches.size();
  ninliers = inliers.size();
  noutliers = nmatches - ninliers;
  ratio = 100.0 * ((float)ninliers / (float)nmatches);

  // Show matching statistics
  std::cout << "Number of Keypoints Image 1: " << nkpts1 << std::endl;
  std::cout << "Number of Keypoints Image 2: " << nkpts2 << std::endl;
  std::cout << "A-KAZE Features Extraction Time (ms): " << takaze << std::endl;
  std::cout << "Matching Descriptors Time (ms): " << tmatch << std::endl;
  std::cout << "Number of Matches: " << nmatches << std::endl;
  std::cout << "Number of Inliers: " << ninliers << std::endl;
  std::cout << "Number of Outliers: " << noutliers << std::endl;
  std::cout << "Inliers Ratio: " << ratio << std::endl << std::endl;

  cimg_library::CImg<float> img1_rgb =
      img1.get_resize(img1.width(), img1.height(), img1.depth(), 3);
  cimg_library::CImg<float> img2_rgb =
      img2.get_resize(img2.width(), img2.height(), img2.depth(), 3);
  draw_keypoints(img1_rgb, kpts1);
  draw_keypoints(img2_rgb, kpts2);
  cimg_library::CImg<float> matched_image;
  draw_matches(img1_rgb, img2_rgb, kpts1, kpts2, inliers, matched_image);
  matched_image.save("../output/matched_images.jpg");
}

/* ************************************************************************* */
int parse_input_options(libAKAZE::AKAZEOptions& options, std::string& img_path1,
                        std::string& img_path2, std::string& homography_path,
                        int argc, char* argv[]) {

  // If there is only one argument return
  if (argc == 1) {
    libAKAZE::show_input_options_help(1);
    return -1;
  }
  // Set the options from the command line
  else if (argc >= 2) {

    // Load the default options
    options = libAKAZE::AKAZEOptions();

    if (!strcmp(argv[1], "--help")) {
      libAKAZE::show_input_options_help(1);
      return -1;
    }

    img_path1 = argv[1];
    img_path2 = argv[2];

    if (argc >= 4) homography_path = argv[3];

    for (int i = 3; i < argc; i++) {
      if (!strcmp(argv[i], "--num_threads")) {
        i = i + 1;
        options.num_threads = atoi(argv[i]);
      } else if (!strcmp(argv[i], "--soffset")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.soffset = atof(argv[i]);
        }
      } else if (!strcmp(argv[i], "--omax")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.omax = atof(argv[i]);
        }
      } else if (!strcmp(argv[i], "--dthreshold")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.dthreshold = atof(argv[i]);
        }
      } else if (!strcmp(argv[i], "--sderivatives")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.sderivatives = atof(argv[i]);
        }
      } else if (!strcmp(argv[i], "--nsublevels")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.nsublevels = atoi(argv[i]);
        }
      } else if (!strcmp(argv[i], "--diffusivity")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.diffusivity = libAKAZE::DIFFUSIVITY_TYPE(atoi(argv[i]));
        }
      } else if (!strcmp(argv[i], "--descriptor")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.descriptor = libAKAZE::DESCRIPTOR_TYPE(atoi(argv[i]));

          if (options.descriptor < 0 || options.descriptor > libAKAZE::MLDB) {
            options.descriptor = libAKAZE::MLDB;
          }
        }
      } else if (!strcmp(argv[i], "--descriptor_channels")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.descriptor_channels = atoi(argv[i]);

          if (options.descriptor_channels <= 0 ||
              options.descriptor_channels > 3) {
            options.descriptor_channels = 3;
          }
        }
      } else if (!strcmp(argv[i], "--descriptor_size")) {
        i = i+1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        }
        else {
          options.descriptor_size = atoi(argv[i]);

          if (options.descriptor_size < 0) {
            options.descriptor_size = 0;
          }
        }
      } else if (!strcmp(argv[i], "--verbose")) {
        options.verbosity = true;
      } else if (!strncmp(argv[i], "--", 2))
        std::cerr << "Unknown command " << argv[i] << std::endl;
    }
  } else {
    std::cerr << "Error introducing input options!!" << std::endl;
    libAKAZE::show_input_options_help(1);
    return -1;
  }

  return 0;
}
