/**
 * @file AKAZE.h
 * @brief Main class for detecting and computing binary descriptors in an
 * accelerated nonlinear scale space
 * @date Oct 07, 2014
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#ifndef AKAZE_SRC_AKAZE_H_
#define AKAZE_SRC_AKAZE_H_

/* ************************************************************************* */
#include <Eigen/Core>

#include "AKAZEConfig.h"
#include "fed.h"
#include "nldiffusion_functions.h"
#include <stdint.h>
/* ************************************************************************* */
namespace libAKAZE {

// Keypoint struct intended to mimic OpenCV.
struct AKAZEKeypoint {
  Eigen::Vector2f pt;
  float size;
  float angle;
  float response;
  int octave;
  int class_id;
};

// This convenience typdef is used to hold binary descriptors such that each bit
// is a value in the descriptor.
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> BinaryVectorX;

struct AKAZEDescriptors {
  std::vector<Eigen::VectorXf> float_descriptor;
  std::vector<BinaryVectorX> binary_descriptor;
};

class AKAZE {
 private:
  // Configuration options for AKAZE
  AKAZEOptions options_;
  // Vector of nonlinear diffusion evolution
  std::vector<TEvolution> evolution_;

  // FED parameters
  // Number of cycles
  int ncycles_;
  // Flag for reordering time steps
  bool reordering_;
  // Vector of FED dynamic time steps
  std::vector<std::vector<float> > tsteps_;
  // Vector of number of steps per cycle
  std::vector<int> nsteps_;

  // Matrices for the M-LDB descriptor computation
  Eigen::MatrixXi descriptorSamples_;
  Eigen::MatrixXi descriptorBits_;

  // Computation times variables in ms
  AKAZETiming timing_;

 public:
  // AKAZE constructor with input options
  // @param options AKAZE configuration options
  // @note This constructor allocates memory for the nonlinear scale space
  AKAZE(const AKAZEOptions& options);

  // Destructor
  ~AKAZE();

  // Allocate the memory for the nonlinear scale space
  void Allocate_Memory_Evolution();

  // This method creates the nonlinear scale space for a given image
  // @param img Input image for which the nonlinear scale space needs to be
  // created
  // @return 0 if the nonlinear scale space was created successfully, -1
  // otherwise
  int Create_Nonlinear_Scale_Space(const RowMatrixXf& img);

  // @brief This method selects interesting keypoints through the nonlinear
  // scale space
  // @param kpts Vector of detected keypoints
  void Feature_Detection(std::vector<AKAZEKeypoint>& kpts);

  // This method computes the feature detector response for the nonlinear scale
  // space
  // @note We use the Hessian determinant as the feature detector response
  void Compute_Determinant_Hessian_Response();

  // This method computes the multiscale derivatives for the nonlinear scale
  // space
  void Compute_Multiscale_Derivatives();

  // This method finds extrema in the nonlinear scale space
  void Find_Scale_Space_Extrema(std::vector<AKAZEKeypoint>& kpts);

  // This method performs subpixel refinement of the detected keypoints fitting
  // a quadratic
  void Do_Subpixel_Refinement(std::vector<AKAZEKeypoint>& kpts);

  // Feature description methods.
  void Compute_Descriptors(std::vector<AKAZEKeypoint>& kpts,
                           AKAZEDescriptors& desc);

  // This method computes the main orientation for a given keypoint
  // @param kpt Input keypoint
  // @note The orientation is computed using a similar approach as described in
  // the original SURF method.
  // See Bay et al., Speeded Up Robust Features, ECCV 2006.
  // A-KAZE uses first order derivatives computed from the nonlinear scale
  // space in contrast to Haar wavelets
  void Compute_Main_Orientation(AKAZEKeypoint& kpt) const;

  // Compute the upright descriptor (not rotation invariant) for the provided
  // keypoint using a
  // rectangular grid similar as the one used in SURF
  // @param kpt Input keypoint
  // @param desc Floating-based descriptor
  // @note Rectangular grid of 20 s x 20 s. Descriptor Length 64. No additional
  // Gaussian weighting is performed. The descriptor is inspired from Bay et
  // al.,
  // Speeded Up Robust Features, ECCV, 2006
  void Get_SURF_Descriptor_Upright_64(const AKAZEKeypoint& kpt,
                                      Eigen::VectorXf& desc) const;

  // Compute the rotation invariant descriptor for the provided keypoint using
  // a
  // rectangular grid similar as the one used in SURF
  // @param kpt Input keypoint
  // @param desc Floating-based descriptor
  // @note Rectangular grid of 20 s x 20 s. Descriptor Length 64. No additional
  // Gaussian weighting is performed. The descriptor is inspired from Bay et
  // al.,
  // Speeded Up Robust Features, ECCV, 2006
  void Get_SURF_Descriptor_64(const AKAZEKeypoint& kpt,
                              Eigen::VectorXf& desc) const;

  // Compute the upright descriptor (not rotation invariant) for the provided
  // keypoint using a
  // rectangular grid similar as the one used in M-SURF
  // @param kpt Input keypoint
  // @param desc Floating-based descriptor
  // @note Rectangular grid of 24 s x 24 s. Descriptor Length 64. The
  // descriptor is inspired
  // from Agrawal et al., CenSurE: Center Surround Extremas for Realtime
  // Feature Detection and Matching,
  // ECCV 2008
  void Get_MSURF_Upright_Descriptor_64(const AKAZEKeypoint& kpt,
                                       Eigen::VectorXf& desc) const;

  // Compute the rotation invariant descriptor for the provided keypoint using
  // a rectangular grid similar as the one used in M-SURF
  // @param kpt Input keypoint
  // @param desc Floating-based descriptor
  // @note Rectangular grid of 24 s x 24 s. Descriptor Length 64. The
  // descriptor is inspired
  // from Agrawal et al., CenSurE: Center Surround Extremas for Realtime
  // Feature Detection and Matching,
  // ECCV 2008
  void Get_MSURF_Descriptor_64(const AKAZEKeypoint& kpt,
                               Eigen::VectorXf& desc) const;

  // Compute the upright (not rotation invariant) M-LDB binary descriptor
  // (maximum descriptor length)
  // @param kpt Input keypoint
  // @param desc Binary-based descriptor
  void Get_Upright_MLDB_Full_Descriptor(const AKAZEKeypoint& kpt,
                                        unsigned char* desc) const;

  // Computes the rotation invariant M-LDB binary descriptor (maximum
  // descriptor length)
  // @param kpt Input keypoint
  // @param desc Binary-based descriptor
  void Get_MLDB_Full_Descriptor(const AKAZEKeypoint& kpt,
                                unsigned char* desc) const;
  /// Compute the upright (not rotation invariant) M-LDB binary descriptor
  /// (specified descriptor length)
  /// @param kpt Input keypoint
  /// @param desc Binary-based descriptor
  void Get_Upright_MLDB_Descriptor_Subset(const AKAZEKeypoint& kpt,
                                          unsigned char* desc);

  /// Computes the rotation invariant M-LDB binary descriptor (specified
  /// descriptor length)
  /// @param kpt Input keypoint
  /// @param desc Binary-based descriptor
  void Get_MLDB_Descriptor_Subset(const AKAZEKeypoint& kpt,
                                  unsigned char* desc);

  // Fill the comparison values for the MLDB rotation invariant descriptor
  void MLDB_Fill_Values(float* values, int sample_step, int level, float xf,
                        float yf, float co, float si, float scale) const;

  // Fill the comparison values for the MLDB upright descriptor
  void MLDB_Fill_Upright_Values(float* values, int sample_step, int level,
                                float xf, float yf, float scale) const;

  // Do the binary comparisons to obtain the descriptor
  void MLDB_Binary_Comparisons(float* values, unsigned char* desc, int count,
                               int& dpos) const;

  // This method saves the scale space into jpg images
  void Save_Scale_Space();

  // This method saves the feature detector responses of the nonlinear scale
  // space into jpg images
  void Save_Detector_Responses();

  // Display timing information
  void Show_Computation_Times() const;

  // Return the computation times
  AKAZETiming Get_Computation_Times() const { return timing_; }
};

/* ************************************************************************* */

// This function sets default parameters for the A-KAZE detector
void setDefaultAKAZEOptions(AKAZEOptions& options);

/// This function computes a (quasi-random) list of bits to be taken
/// from the full descriptor. To speed the extraction, the function creates
/// a list of the samples that are involved in generating at least a bit
/// (sampleList)
/// and a list of the comparisons between those samples (comparisons)
/// @param sampleList
/// @param comparisons The matrix with the binary comparisons
/// @param nbits The number of bits of the descriptor
/// @param pattern_size The pattern size for the binary descriptor
/// @param nchannels Number of channels to consider in the descriptor (1-3)
/// @note The function keeps the 18 bits (3-channels by 6 comparisons) of the
/// coarser grid, since it provides the most robust estimations
void generateDescriptorSubsample(Eigen::MatrixXi& sampleList,
                                 Eigen::MatrixXi& comparisons, int nbits,
                                 int pattern_size, int nchannels);

// This function checks descriptor limits for a given keypoint
inline void check_descriptor_limits(int& x, int& y, int width, int height);

// This function computes the value of a 2D Gaussian function
inline float gaussian(float x, float y, float sigma) {
  return expf(-(x * x + y * y) / (2.0f * sigma * sigma));
}

// This funtion rounds float to nearest integer
inline int fRound(float flt) { return (int)(flt + 0.5f); }

}  // namespace libAKAZE

#endif  // AKAZE_SRC_AKAZE_H_
