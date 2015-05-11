.. highlight:: c++

.. default-domain:: cpp

.. _documentation-features:

========
Features
========

Feature detection and description is a major area of focus in Computer Vision. While SIFT remains the gold standard because of its robustness and matching performance, many other detectors and descriptors are used and often have other competitive advantages. Theia presents friendly classes for feature detection and decription such that the interface is always the same regardless of the methods used. Note that all keypoint and descriptor extraction methods we perform automatic conversion to grayscale images if necessary.

.. NOTE:: The keypoint detection and descriptor extraction methods here are generally not thread-safe. Use with caution when multi-threading.

Keypoints
=========

.. class:: Keypoint

The base :class:`Keypoint` class is a glorified struct that holds information about a keypoint that has been detected with a :class:`KeypointDetector`. Information about the keypoint's position, strength, scale, and orientation can be easily added and retrieved. The type of keypoint can be retrieved with the :func:`keypoint_type()` function.

	.. code-block:: c++

           class Keypoint {
	    public:
	      enum KeypointType {
	          INVALID = -1,
		  OTHER = 0,
		  SIFT,
		  AGAST,
		  BRISK
		  };

	      Keypoint(double x, double y, KeypointType type);
	      ~Keypoint() {}

	      // Required Keypoint type.
	      inline KeypointType keypoint_type() const;
	      inline void set_keypoint_type(KeypointType type);

	      // Required Variable x.
	      inline double x() const;
	      inline void set_x(double x);

	      // Required Variable y.
	      inline double y() const;
	      inline void set_y(double y);

	      // Optional variable strength.
	      inline bool has_strength() const;
	      inline double strength() const;
	      inline void set_strength(double strength);

	      // Optional variable scale.
	      inline bool has_scale() const;
	      inline double scale() const;
	      inline void set_scale(double scale);

	      // Optional variable orientation.
	      inline bool has_orientation() const;
	      inline double orientation() const;
	      inline void set_orientation(double orientation);
	   };

Keypoint Detector
=================

.. class:: KeypointDetector

Detecting keypoints with Theia is very simple, and we have implemented a number of keypoint detectors that are commonly used in Computer Vision. Each keypoint detector is derived from the virtual class :class:`KeypointDetector`. Each derived class must implement the :func:`DetectKeypoints` method

.. function:: bool KeypointDetector::Initialize()

    This method initializes any internal parameters that must be generated,
    precalculated, or otherwise are independent of the image. The
    :func:`Initialize()` function must be called before using the keypoint
    detector.

.. function:: bool KeypointDetector::DetectKeypoints(const FloatImage& input_image, std::vector<Keypoint>* output_keypoints)

  ``input_image``: The image that you want to detect keypoints on.

  ``ouput_keypoints``: A pointer to a vector that will hold the keypoints
    detected. Note that the vector should be empty when passed to this
    function. The caller is responsible for deleting the keypoints.

  .. code-block:: c++

    // Assume var keypoint_detector was created with one of the constructors below.

    FloatImage input_image(input_image_filename);
    const bool initialization_success = keypoint_detector.Initialize();

    // Container for the detected keypoints.
    std::vector<Keypoint> detected_keypoint;
    const bool detection_success =
        keypoint_detector.DetectKeypoints(input_image, &detected_keypoints);


The following keypoint detectors have been implemented in Theia (class constructors are given):

.. function:: SiftDetector::SiftDetector(int num_octaves, int num_scale_levels, int first_octave)

    The algorithm originally proposed by [Lowe]_ that uses the `VLFeat
    <http://www.vlfeat.org>`_ as the underlying engine.

    Specify the number of image octaves, number of scale levels per octave, and
    where the first octave should start. The default constructor sets these values
    to values -1 (i.e., as many octaves as can be generated), 3, and 0 (i.e., the
    source image)

.. function:: AgastDetector::AgastDetector(AstPattern pattern, int threshold, bool nonmax_suppression)

    The improved FAST detection scheme of [Mair]_ et al.

    ``enum AstPattern`` specifies one of 4 types of sampling patterns for the
    AGAST corner detect: ``AGAST5_8`` is the AGAST pattern with an 8 pixel mask,
    ``AGAST7_12D`` is the AGAST diamond pattern with a 12 pixel mask,
    ``AGAST7_12S`` is the square configuration, and ``OAST9_16`` is the 16 pixel
    mask. By default, we the detector uses ``AGAST5_8`` with a threshold of 30 and
    nonmaximum suppression turn on. More details on the configurations can be
    found at the `AGAST Project website
    <http://www6.in.tum.de/Main/ResearchAgast>`_

.. function:: BriskDetector::BriskDetector(int threshold, int num_octaves)

  The "Binary Robust Invariant Scalable Keypoints" algorithm of [Leutenegger]_
  et al.

  Specify the threshold for keypoint scores (default is 30) and the number of
  octaves to downsample the image (default is 3).

Descriptors
===========

Theia uses a semi-generic interface for all descriptor types, namely, floating point and binary descriptors. For floating point descriptors (e.g., SIFT) we use Eigen::VectorXf and set the number of entries to equal the dimension of the descriptor. This way, we can utilize Eigen's speed and optimizations to get the most efficient and accurate representation of the descriptors. For binary descriptors, we define a new type in the Eigen namespace: ``Eigen::BinaryVectorX``. This vector is a custom type (defined in theia/alignment/alignment.h) that holds binary descriptors such that each bit corresponds to the descriptor dimension. This allows for the same interface between float and binary descriptors, while still utilizing the efficiency of SSE instructions when available.

DescriptorExtractor
===================

.. class:: DescriptorExtractor

  We enforce a :class:`DescriptorExtractor` interface similar to the
  :class:`KeypointDetector` so that we can extract descriptors at runtime. Each
  descriptor has a corresponding extractor class that is used to compute that
  descriptor given keypoints. However, we must call the :func:`Initialize()`
  method before computing descriptors.

.. function:: bool DescriptorExtractor::Initialize()

  This method initializes any internal parameters that must be generated,
  precalculated, or otherwise are independent of the image. The
  :func:`Initialize()` function must be called before using the descriptor
  extractor.

.. function:: bool DescriptorExtractor::ComputeDescriptor(const FloatImage& input_image, const Keypoint& keypoint, Eigen::VectorXf* float_descriptor)
.. function:: bool DescriptorExtractor::ComputeDescriptor(const FloatImage& input_image, const Keypoint& keypoint, Eigen::BinaryVectorXf* binary_descriptor)

  This method computes the descriptor of a single keypoint.

  ``input_image``: The image that you want to detect keypoints on.

  ``keypoint``: The keypoint that the descriptor will be computed from.

  ``float_descriptor or binary_descriptor``: The descriptor computed for the
  given keypoint.

  ``returns``: True on if the descriptor was extracted, false otherwise.

.. function:: bool DescriptorExtractor::ComputeDescriptors(const FloatImage& input_image, std::vector<Keypoint>* keypoints, std::vector<Eigen::VectorXf>* float_descriptors)
.. function:: bool DescriptorExtractor::ComputeDescriptors(const FloatImage& input_image, std::vector<Keypoint>* keypoints, std::vector<Eigen::BinaryVectorXf>* binary_descriptors)

    Compute many descriptors from the input keypoints. Note that not all
    keypoints are guaranteed to result in a descriptor. Only valid descriptors
    (and feature positions) are returned in the output parameters.

    ``input_image``: The image that you want to detect keypoints on.

    ``keypoints``: An input vector of the keypoint pointers that will have
    descriptors extracted. Keypoints that were not able to have a descriptor
    extracted are removed.

    ``float_descriptors or binary_descriptors``: A container for the descriptors
    that have been created based on the type of descriptor that is being
    extracted. Eigen::VectorXf is used for extracting float descriptors (e.g.,
    SIFT) while Eigen::BinaryVectorX is used for float descriptors.

.. function:: bool DescriptorExtractor::DetectAndExtractDescriptors(const FloatImage& input_image, std::vector<Keypoint>* keypoints, std::vector<Eigen::VectorXf>* float_descriptors)
.. function:: bool DescriptorExtractor::DetectAndExtractDescriptors(const FloatImage& input_image, std::vector<Keypoint>* keypoints, std::vector<Eigen::BinaryVectorXf>* binary_descriptors)

    Detects keypoints and extracts descriptors using the default keypoint
    detector for the corresponding descriptor. For SIFT, this is the SIFT
    keypoint detector, and for BRIEF, BRISK, and FREAK this is the BRISK
    keypoint detector. This has the potential to be faster because it may avoid
    recomputing certain member variables.

    ``input_image``: The image that you want to detect keypoints on.

    ``keypoints``: An output vector of the keypoint points that have been
    detected and successfully had descriptors extracted.

    ``float_descriptors or binary_descriptors``: A container for the descriptors
    that have been created based on the type of descriptor that is being
    extracted. Eigen::VectorXf is used for extracting float descriptors (e.g.,
    SIFT) while Eigen::BinaryVectorX is used for float descriptors.

  .. code-block:: c++

    // Open image we want to extract features from.
    FloatImage input_image(input_image_filename);

    // Detect keypoints.
    SiftDetector sift_keypoint_detector;
    bool keypoint_detector_init = sift_keypoint_detector.Initialize();
    const bool keypoint_init_success = sift_keypoint_detector.Initialize();
    std::vector<Keypoint> sift_keypoints;
    const bool detection_success =
        sift_keypoint_detector.DetectKeypoints(input_image, &sift_keypoints);

    // Initialize descriptor extractor.
    SiftDescriptorExtractor sift_extractor;
    const bool descriptor_init_succes = sift_extractor.Initialize();

    // E.g., compute a single descriptor
    Eigen::VectorXf sift_descriptor;
    bool sift_success =
      sift_extractor.ComputeDescriptor(input_image, keypoint[0], &sift_descriptor);

    // E.g., compute many descriptors.
    std::vector<Eigen::VectorXf> sift_descriptors;
    const bool extraction_success =
      sift_extractor.ComputeDescriptors(image, &sift_keypoints, &sift_descriptors)

We implement the following descriptor extractors (and corresponding descriptors)
in Theia (constructors are given).

.. class:: SiftDescriptorExtractor

.. function:: SiftDescriptorExtractor::SiftDescriptorExtractor(int num_octaves, int num_scale_levels, int first_octave)

  The algorithm originally proposed by [Lowe]_ that uses the `VLFeat
  <http://www.vlfeat.org>`_ as the underlying engine.

  We only implement the standard 128-dimension descriptor. Specify the number
  of image octaves, number of scale levels per octave, and where the first
  octave should start. The default constructor sets these values to values -1
  (i.e., as many octaves as can be generated), 3, and 0 (i.e., the source
  image). Typically these parameters are set to match the :class:`SiftDetector`
  parameters.

.. NOTE:: This algorithm is patented and commercial use requires a license.

.. class:: BriefDescriptorExtractor

.. function:: BriefDescriptorExtractor::BriefDescriptorExtractor(int patch_sample_size, const int num_bytes)

   The [BRIEF]_ algorithm is a binary algorithm that operates on local image
   patches around a keypoint or a point of interest. The binary values are set
   by randomly choosing two pixels to compare within the patch. The same random
   pattern must be used in order to compare BRIEF descriptors (each
   :class:`BriefDescriptorExtractor` object creates exactly one pattern that
   may be used repeatedly).


.. class:: FreakDescriptorExtractor

.. NOTE:: This algorithm is currently unstable. Further testing is required.

.. function:: FreakDescriptorExtractor::FreakDescriptorExtractor(bool rotation_invariant, bool scale_invariant, int num_octaves)

  The "Fast Retina Keypoint" algorithm for binary descriptors proposed by [Alahi]_ et al.

  ``rotation_invariant``: Set to true if you want to normalize the orientation of the keypoints before computing the descriptor.

  ``scale_invariant``: Set to true if you want to normalize the scale of keypoints before computing the descriptor.

  ``num_octaves``: The number of octaves that the keypoints span.

  The :class:`FreakDescriptorExtractor` is typically used with the
  :class:`BriskDetector` to detect keypoints.

.. class:: BriskDescriptorExtractor

.. NOTE:: This algorithm is currently unstable. Further testing is required.

.. function:: BriskDescriptorExtractor::BriskDescriptorExtractor(bool rotation_invariant, bool scale_invariant, float pattern_scale)

  The "Binary Robust Invariant Scalable Keypoints" algorithm for binary descriptors of [Leutenegger]_
  et al.

  ``rotation_invariant``: Set to true if you want to normalize the orientation of the keypoints before computing the descriptor.

  ``scale_invariant``: Set to true if you want to normalize the scale of keypoints before computing the descriptor.

  ``pattern_scale``: Scale of the BRISK pattern to use.


Feature Matching
================

Features are useful in SfM because they can provide sparse matches between
images, which can then provide geometric constrainst for the poses between these
images. As such, feature matching is a very critical process in the context of
multi-view geometry. We provide a generic interface for feature matching that
works with binary descriptors or float descriptors.

.. class:: FeatureMatcher

The :class:`FeatureMatcher` is templated on a :class:`DistanceMetric` that
describes how to compute the distance between two matches (we provide L2 and
Hamming). The matcher is intended for all-to-all matching for SfM reconstruction.

.. function:: void FeatureMatcher::AddImage(const std::vector<Keypoint>* keypoints, const std::vector<DescriptorType>* descriptors)

  Adds an image to the matcher with no known intrinsics for this image. The
  caller still owns the keypoints and descriptors so they must remain valid
  objects throughout the matching.

.. function:: void FeatureMatcherAddImage(const std::vector<Keypoint>* keypoints, const std::vector<DescriptorType>* descriptors, const CameraIntrinsics& intrinsics)

  Adds an image to the matcher with the known camera intrinsics. The
  intrinsics (if known) are useful for geometric verification. The caller
  still owns the keypoints and descriptors so they must remain valid objects
  throughout the matching.

.. function:: void FeatureMatcher::MatchImages(const FeatureMatcherOptions& matcher_options, std::vector<ImagePairMatch>* matches)

  Matches features between all images. No geometric verification is
  performed. Only the matches which pass the have greater than
  min_num_feature_matches are returned.

.. function:: void FeatureMatcher::MatchImagesWithGeometricVerification(const FeatureMatcherOptions& matcher_options, const VerifyTwoViewMatchesOptions& verification_options, std::vector<ImagePairMatch>* matches)

  Matches features between all images. Only the matches that pass the
  geometric verification are returned. Camera intrinsics are used for
  geometric verification if the image was added with known intrinsics.

.. NOTE:: This method is tuned specifically for image to image matching and is only
   applicable to float descriptors such as SIFT.

.. class:: ImagePairMatch

Matches are defined as feature coordinates between two image. If geometric
verification is performed then the two-view geometry is also specified and the
returned matches are only the inlier matches after geometric verification.

.. member:: int ImagePairMatch::image1_index
.. member:: int ImagePairMatch::image2_index

  The index of the current image pair that has been matched. This index is
  relative to the order that images were input with the
  :func:`FeatureMatcher::AddImage` method.

.. member:: TwoViewInfo FeatureMatcher::twoview_info

  If geometric verification is performed, then the ``twoview_info`` describes
  the two-view geometry (i.e., relative pose) between the two images.

.. member:: std::vector<:class:`FeatureCorrespondence`> FeatureMatcher::correspondences

  A :class:`FeatureCorrespondence` contains two Eigen::Vector2d's named
  feature1, and feature2. These represent the image coordinates of the matched
  features. If geometric verification is performed then these features are the
  inlier features.

.. class:: FeatureMatcherOptions

  The options specified for feature matching. Adjusting these optiosn will
  change the number of matched features as well as the quality for matching.

.. member:: int FeatureMatcherOptions::num_threads

  DEFAULT: ``1``

  The number of threads to use for image-to-image matching. The more threads
  used, the faster the matching will be.

.. member:: bool FeatureMatcherOptions::keep_only_symmetric_matches

  DEFAULT: ``true``

  The quality of feature matching can be greatly improved by only keeping
  matches that are mutual. That is, for feature ``x`` in image 1 and feature
  ``y`` in image 2, a high quality match is formed when ``y`` is the best match
  for ``x`` and ``x`` is also the best match for ``y``. When
  ``keep_only_symmetric_matches`` is enabled, only mutual matches are considered
  valid.

.. member:: bool FeatureMatcherOptions::use_lowes_ratio

  DEFAULT: ``true``

.. member:: float FeatureMatcherOptions::lowes_ratio

  DEFAULT: ``0.8``

  Good feature matches should be very apparent. That is, the best match for a
  given feature should be much better than all other candidate matches for a
  given feature. Lowes ratio is defined as the ratio between the top match
  distance and the second best match distance. If this ratio is higher than
  ``lowes_ratio`` then that means that the top match is not much better than the
  second best match. If ``use_lowes_ratio`` is set to ``true`` then only the
  feature matches which pass the Lowes ratio test are kept.

.. member:: int FeatureMatcherOptions::min_num_feature_matches

  DEFAULT: ``30``

  Images are only considered to be successfully matched if they contain a
  sufficient number of feature matches between them. ``min_num_feature_matches``
  is the minimum number of valid feature matches (or verified matches) that must
  exist between two images in order to consider the matches as valid. All other
  matches are considered failed matches and are not added to the output.

Matching Strategies
-------------------

We have implemented two types of :class:`FeatureMatcher` with the interface described above.

.. class:: BruteForceFeatureMatcher

Matches are computed using an exhausitve brute force search through all
matches. The search is the slowest but has the highest accuracy.

.. class:: CascadeHashingFeatureMatcher

Features are matched through a cascade hashing approach as described by
[Cheng]_. Hash tables with extremely fast lookups are created without needing to
train the data, resulting in an extremely fast and accurate matcher. This is the
recommended approach for matching image sets.


Using the feature matcher
-------------------------

The intended use for these classes is for matching photos in image collections,
so all pairwise matches are computed. Matching with geometric verification is
also possible. Typical use case is:

.. code-block:: c++

      FeatureMatcher matcher;
      for (int i = 0; i < num_images_to_match; i++) {
        matcher.AddImage(keypoints[i], descriptors[i]);

       // Or, you could add the image with known intrinsics for use during
       // geometric verification.
        matcher.AddImage(keypoints[i], descriptors[i], intrinsics[i]);
      }
      std::vector<ImagePairMatch> matches;
      FeatureMatcherOptions matcher_options;
      matcher.MatchImages(matcher_options, &matches);
          Or, with geometric verification:
      VerifyTwoViewMatchesOptions geometric_verification_options;
      matcher.MatchImages(match_options,
                          geometric_verification_options,
                          &matches);
