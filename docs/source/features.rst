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

Descriptors
===========

Theia uses a semi-generic interface for all descriptor types. For floating point descriptors (e.g., SIFT) we use Eigen::VectorXf and set the number of entries to equal the dimension of the descriptor. This way, we can utilize Eigen's speed and optimizations to get the most efficient and accurate representation of the descriptors.

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

  This method computes the descriptor of a single keypoint.

  ``input_image``: The image that you want to detect keypoints on.

  ``keypoint``: The keypoint that the descriptor will be computed from.

  ``float_descriptor``: The descriptor computed for the
  given keypoint.

  ``returns``: True on if the descriptor was extracted, false otherwise.

.. function:: bool DescriptorExtractor::ComputeDescriptors(const FloatImage& input_image, std::vector<Keypoint>* keypoints, std::vector<Eigen::VectorXf>* float_descriptors)

    Compute many descriptors from the input keypoints. Note that not all
    keypoints are guaranteed to result in a descriptor. Only valid descriptors
    (and feature positions) are returned in the output parameters.

    ``input_image``: The image that you want to detect keypoints on.

    ``keypoints``: An input vector of the keypoint pointers that will have
    descriptors extracted. Keypoints that were not able to have a descriptor
    extracted are removed.

    ``float_descriptors``: A container for the descriptors
    that have been created based on the type of descriptor that is being
    extracted.

.. function:: bool DescriptorExtractor::DetectAndExtractDescriptors(const FloatImage& input_image, std::vector<Keypoint>* keypoints, std::vector<Eigen::VectorXf>* float_descriptors)

    Detects keypoints and extracts descriptors using the default keypoint
    detector for the corresponding descriptor. For SIFT, this is the SIFT
    keypoint detector. This has the potential to be faster because it may avoid
    recomputing certain member variables.

    ``input_image``: The image that you want to detect keypoints on.

    ``keypoints``: An output vector of the keypoint points that have been
    detected and successfully had descriptors extracted.

    ``float_descriptors``: A container for the descriptors
    that have been created based on the type of descriptor that is being
    extracted. Eigen::VectorXf is used for extracting float descriptors (e.g.,
    SIFT).

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


Feature Matching
================

Features are useful in SfM because they can provide sparse matches between
images, which can then provide geometric constrainst for the poses between these
images. As such, feature matching is a very critical process in the context of
multi-view geometry. We provide a generic interface for feature matching that
works with binary descriptors or float descriptors.

For feature matching, we implement an abstract :class:`FeatureMatcher` class that
serves as an abstract class for various feature-matching methods. The
:class:`FeatureMatcher` class takes keypoints, descriptors, and optionally
camera intrinsics (if known) and performs all-pairs feature matching between images.

.. class:: FeatureMatcher

The :class:`FeatureMatcher` is templated on a :class:`DistanceMetric` that
describes how to compute the distance between two matches (we provide L2 and
Hamming). The matcher is intended for all-pairs image matching for SfM
reconstruction.

.. function:: FeatureMatcher::FeatureMatcher(const FeatureMatcherOptions& options)

   Initializes a feature matcher based on the options.

.. function:: void FeatureMatcher::AddImage(const std::string& image_name, const std::vector<Keypoint>& keypoints, const std::vector<DescriptorType>& descriptors)

  Adds an image to the matcher with no known intrinsics for this image. The
  image name must be a unique identifier.

.. function:: void FeatureMatcherAddImage(const std::string& image_name, const std::vector<Keypoint>& keypoints, const std::vector<DescriptorType>& descriptors, const CameraIntrinsics& intrinsics)

  Adds an image to the matcher with the known camera intrinsics. The intrinsics
  (if known) are used for geometric verification. The image name must be a
  unique identifier.

.. function:: void FeatureMatcher::MatchImages(std::vector<ImagePairMatch>* matches)

  Matches features between all images. No geometric verification is
  performed. Only successful image matches will be returned.

.. function:: void FeatureMatcher::MatchImagesWithGeometricVerification(const VerifyTwoViewMatchesOptions& verification_options, std::vector<ImagePairMatch>* matches)

  Matches features between all images. Only the matches that pass the
  geometric verification are returned. Camera intrinsics are used for
  geometric verification if the image was added with known intrinsics.

.. function:: void FeatureMatcher::SetImagePairsToMatch(const std::vector<std::pair<std::string, std::string> >& pairs_to_match)

  Set the image pairs that will be matched when MatchImages or
  MatchImagesWithGeometricVerification is called. This is an optional method; if
  it is not called, then all possible image-to-image pairs will be matched. The
  vector should contain unique pairs of image names that should be matched.


Feature Matching Options
------------------------

Theia allows for a variety of parameters to be tuned for feature
matching. Setting these parameters will have an effect on things such as
matching performance, efficiency, memory, and more.

.. class:: FeatureMatcherOptions

  The options specified for feature matching. Adjusting these optiosn will
  change the number of matched features as well as the quality for matching.

.. member:: int FeatureMatcherOptions::num_threads

  DEFAULT: ``1``

  The number of threads to use for image-to-image matching. The more threads
  used, the faster the matching will be.

.. member:: bool FeatureMatcherOptions::match_out_of_core

  DEFAULT: ``false``

  Matching can be performed out-of-core or all in memory. For large datasets, it
  is advisable to utilize the out-of-core matching. This strategy will save
  features to disk and utilize an LRU cache to minimize disk IO and take
  advantage of cache-locality.

.. member:: std::string FeatureMatcherOptions::keypoints_and_descriptors_output_dir

  DEFAULT: ``""``

  If out-of-core matching is enabled, this is the directory where features will
  be written to and read from disk.

.. member:: int FeatureMatcherOptions::cache_capacity

  DEFAULT: ``128``

  If out-of-core matching is enabled, this is the maximum number of images to
  store in the cache at a given time. The larger this number, the more memory is
  required for matching.

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


Output of Feature Matching
--------------------------

 The output of the matching process is a vector of :class:`ImagePairMatch`. Each
 :class:`ImagePairMatch` contains matching information for feature matches
 between two views.

.. class:: ImagePairMatch

Matches are defined as feature coordinates between two image. If geometric
verification is performed then the two-view geometry is also specified and the
returned matches are only the inlier matches after geometric verification.

.. member:: std::string ImagePairMatch::image1
.. member:: std::string ImagePairMatch::image2

  The unique names of the current image pair that have been matched.

.. member:: TwoViewInfo ImagePairMatch::twoview_info

  If geometric verification is performed, then the ``twoview_info`` describes
  the two-view geometry (i.e., relative pose) between the two images.

.. member:: std::vector<FeatureCorrespondence> ImagePairMatch::correspondences

  A :class:`FeatureCorrespondence` contains two feature locations named
  feature1, and feature2. These represent the image coordinates of the matched
  features. If geometric verification is performed then these features are the
  inlier features.


Using the feature matcher
-------------------------

We have implemented two types of :class:`FeatureMatcher` with the interface described above.

.. class:: BruteForceFeatureMatcher

  Matches are computed using an exhausitve brute force search through all
  matches. The search is the slowest but has the highest accuracy.

.. class:: CascadeHashingFeatureMatcher

  Features are matched through a cascade hashing approach as described by
  [Cheng]_. Hash tables with extremely fast lookups are created without needing to
  train the data, resulting in an extremely fast and accurate matcher. This is the
  recommended approach for matching image sets.


The intended use for the :class:`FeatureMatcher` is for matching photos in image collections,
so all pairwise matches are computed. Typical use case is:


.. code-block:: c++

      FeatureMatcherOptions matcher_options;
      BruteForceFeatureMatcher matcher(matcher_options);
      // Or to instantiate the cascade hashing matcher:
      CascadeHashingFeatureMatcher matcher(matcher_options);

      // Add image features to the matcher.
      for (int i = 0; i < num_images_to_match; i++) {
        matcher.AddImage(image_name[i], keypoints[i], descriptors[i]);

       // Or, you could add the image with known intrinsics for use during
       // geometric verification.
        matcher.AddImage(image_name[i], keypoints[i], descriptors[i], intrinsics[i]);
      }
      std::vector<ImagePairMatch> matches;
      matcher.MatchImages(&matches);

      // Or, with geometric verification:
      VerifyTwoViewMatchesOptions geometric_verification_options;
      matcher.MatchImages(geometric_verification_options, &matches);

By adjusting the :class:`FeatureMatcherOptions` (described above) you can
control various setting such as multithreading, in-core vs out-of-core, etc. The
outpute of the matching process is a vector of :class:`ImagePairMatch`. Each
:class:`ImagePairMatch` contains matching information for feature matches
between two views.


Implementing a New Matching Strategy
------------------------------------

Given that the :class:`FeatureMatcher` class is an abstract interface,
implementing a new matching strategy is extremely simple. The simplest way to do
this is to derive a new class from the :class:`FeatureMatcher` class and
implement the protected method :func:`MatchImagePair`

.. function:: bool FeatureMatcher::MatchImagePair(const KeypointsAndDescriptors& features1, const KeypointsAndDescriptors& features2, std::vector<FeatureCorrespondence>* matched_features)

   This protected function takes in two sets of features and outputs the feature
   matches between them. This is a pure virual function in the
   :class:`FeatureMatcher` class and must be implemented by any derived
   classes. For instance, the :class:`BruteForceFeatureMatcher` implements this
   method by computing the pairwise distance between all features and choosing
   the correspondences as the features with the smallest distance between them.

   When implementing this method in a derived class you will automatically get
   all of the great benefits of the abstract :class:`FeatureMatcher` class
   without having to explicitly write code to handle them. These benefits include:

   * Multithreaded matching
   * Ability to utilize out-of-core matching
   * Optional geometric verification

For examples on how to implemente new matchers as derived classes, check out the
:class:`BruteForceFeatureMatcher` implementation.
