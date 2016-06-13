//=============================================================================
//
// akaze_features.cpp
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
 * @file akaze_features.cpp
 * @brief Main program for detecting and computing binary descriptors in an
 * accelerated nonlinear scale space
 * @date Oct 07, 2014
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#include <ctime>

#include "cimg/CImg.h"
#include "src/AKAZE.h"
#include "src/AKAZEConfig.h"
#include "src/utils.h"
#include "timer/timer.hpp"

/* ************************************************************************* */
/**
 * @brief This function parses the command line arguments for setting A-KAZE
 * parameters
 * @param options Structure that contains A-KAZE settings
 * @param img_path Path for the input image
 * @param kpts_path Path for the file where the keypoints where be stored
 */
int parse_input_options(libAKAZE::AKAZEOptions& options, std::string& img_path,
                        std::string& kpts_path, int argc, char* argv[]);

/* ************************************************************************* */
int main(int argc, char* argv[]) {
  // Variables
  libAKAZE::AKAZEOptions options;
  std::string img_path, kpts_path;

  // Variable for computation times.
  double tdet = 0.0, tdesc = 0.0;

  // Parse the input command line options
  if (parse_input_options(options, img_path, kpts_path, argc, argv)) {
    return -1;
  }

  if (options.verbosity) {
    std::cout << "Check AKAZE options:" << std::endl;
    std::cout << options << std::endl;
  }

  // Try to read the image and if necessary convert to grayscale. CImg will
  // throw an error and crash if the image could not be read.
  cimg_library::CImg<float> img(img_path.c_str());
  RowMatrixXf img_32;
  libAKAZE::ConvertCImgToEigen(img, img_32);
  img_32 /= 255.0;

  // Don't forget to specify image dimensions in AKAZE's options.
  options.img_width = img_32.cols();
  options.img_height = img_32.rows();

  // Extract features.
  std::vector<libAKAZE::AKAZEKeypoint> kpts;
  timer::Timer timer;
  libAKAZE::AKAZE evolution(options);
  evolution.Create_Nonlinear_Scale_Space(img_32);
  evolution.Feature_Detection(kpts);
  tdet = timer.elapsedMs();

  // Compute descriptors.
  libAKAZE::AKAZEDescriptors desc;
  timer.reset();
  evolution.Compute_Descriptors(kpts, desc);
  tdesc = timer.elapsedMs();

  // Summarize the computation times.
  evolution.Show_Computation_Times();
  evolution.Save_Scale_Space();
  std::cout << "Number of points: " << kpts.size() << std::endl;
  std::cout << "Time Detector: " << tdet << " ms" << std::endl;
  std::cout << "Time Descriptor: " << tdesc << " ms" << std::endl;

  // Save keypoints in ASCII format.
  if (!kpts_path.empty()) {
    if (options.descriptor < libAKAZE::MLDB_UPRIGHT) {
      libAKAZE::save_keypoints(kpts_path, kpts, desc.float_descriptor, true);
    } else {
      libAKAZE::save_keypoints(kpts_path, kpts, desc.binary_descriptor, true);
    }
  }

  // Convert the input image to RGB.
  cimg_library::CImg<float> rgb_image =
      img.get_resize(img.width(), img.height(), img.depth(), 3);
  libAKAZE::draw_keypoints(rgb_image, kpts);
  rgb_image.save("../output/detected_features.jpg");
}

/* ************************************************************************* */
int parse_input_options(libAKAZE::AKAZEOptions& options, std::string& img_path,
                        std::string& kpts_path, int argc, char* argv[]) {

  // If there is only one argument return
  if (argc < 2) {
    libAKAZE::show_input_options_help(0);
    return -1;
  }
  // Set the options from the command line
  else if (argc >= 3) {
    options = libAKAZE::AKAZEOptions();
    kpts_path = "./keypoints.txt";

    if (!strcmp(argv[1], "--help")) {
      libAKAZE::show_input_options_help(0);
      return -1;
    }

    img_path = argv[1];

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
        } else
          options.nsublevels = atoi(argv[i]);
      } else if (!strcmp(argv[i], "--diffusivity")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else
          options.diffusivity = libAKAZE::DIFFUSIVITY_TYPE(atoi(argv[i]));
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
      } else if (!strcmp(argv[i],"--descriptor_size")) {
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
      } else if (!strcmp(argv[i], "--save_scale_space")) {
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else {
          options.save_scale_space = (bool)atoi(argv[i]);
        }
      } else if (!strcmp(argv[i], "--verbose")) {
        options.verbosity = true;
      } else if (!strcmp(argv[i], "--output")) {
        options.save_keypoints = true;
        i = i + 1;
        if (i >= argc) {
          std::cerr << "Error introducing input options!!" << std::endl;
          return -1;
        } else
          kpts_path = argv[i];
      }
    }
  }

  return 0;
}
