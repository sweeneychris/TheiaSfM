.. _chapter-faq:

============================
FAQs, and General Guidelines
============================

Answers to common questions about SfM and Theia, as well as some guidelines for
obtaining good 3D reconstructions.

Building
========

#. Does Theia compile on Windows?

    Theia does not work out-of-the-box on most Windows machines,
    unfortunately. This is largely due to the fact that Theia uses c++11
    features that some versions of Visual Studio do not implement. However, many
    users have successfully patched Theia to work with Windows without too much
    trouble. Please feel free to submit GitHub issues or patches to document
    which lines of code are causing compiler errors on Windows and I will gladly
    incorporate the fixes into the library.

#. Does Theia work on Android?

    Not currently. This is largely due to the fact that Theia extensively uses
    the Glog (google logging) library, which does not currently work on
    Android. Further, some of the sparse linear solvers that Theia uses cannot
    be used on Android.

Creating 3D Reconstructions
===========================

#. How do I get a good 3D reconstruction?

    Obtaining a good 3D reconstruction is not always trivial, but there are a
    few things that can help the Structure-from-Motion process converge to a better solution.

    #. Use accurate camera calibration (camera intrinsic paramters). Theia will
       attempt to determine camera focal lengths from EXIF, but the best results
       are obtained with explicit calibration. If calibration is accurate then
       you may not need to optimize camera intrinsics during bundle adjustment.

    #. Ensure the images have good parallax. In order to triangulate 3D point
       accurately, images must be sufficiently far apart so that the
       triangulation is well-constrained. Scenes where the cameras are near
       rotation-only motions are usually difficult to reconstruct.

    #. Ensure there are a sufficient number of features extracted and
       matched. Usually extracting several thousand features (for a 1MP image)
       and having two-view matches with several hundred or thousand matches will
       results in good reconstructions.

#. Incremental vs Global SfM

    The difference between incremental and global SfM pipelines may be found in
    the :ref:`chapter-sfm` chapter. Generally speaking, incremental pipelines
    are slower (less scalable), but more robust and accurate. This is due to the
    extensive use of RANSAC for outlier filtering and repeated use of bundle
    adjustment. Global SfM is more scalable and can often converge to a pretty
    good solution, but is more susceptible to outliers. If global SfM is not
    working well for you and reconstruction time is not an issue then
    incremental SfM is recommended.

#. I ran SfM but cannot obtain a good reconstruction. What paramters should I change?

   See the above question about how to get a good reconstruction. Useful
   parameters to tweak are feature thresholds (try to get more features),
   increasing lowes ratio for matching (try to get more matches), increasing the
   maximum reprojection error allowed (more points used for bundle adjustment),
   and increasing the number of retriangulation iterations for global sfm.


#. I have a good 3D reconstruction... now what?

   There are lots of interesting things you can do with 3D reconstructions!  For
   instance, image-based localization, scene recognition/understanding,
   multi-view stereo, model-based tracking, and more. Most of these applications
   are not implemented within Theia, but the algorithms within Theia will help
   you create such applications.

Extending Theia
===============

#. How do I implement a new SfM pipeline?

    New SfM pipelines can be easily created by deriving a new class from the
    :class:`ReconstructionEstimator` class. Most likely, though, you will want
    to modify the existing incremental or global SfM pipelines to suit your
    needs.

#. How do I implement a new rotation or position estimation algorithm for Global SfM?

    Theia utilizes a modular approach to global SfM estimation. Our global SfM
    pipeline calls rotation and position estimation methods through the abstract
    :class:`RotationEstimator` and :class:`PositionEstimator` classes (see
    :ref:`chapter-sfm`). To implement a new rotation or position estimation
    algorithm, simply derive from one of these classes and add your new method
    to either :func:`GlobalReconstructionEstimator::EstimateGlobalRotations()`
    or :func:`GlobalReconstructionEstimator::EstimatePositions()`. Then you can
    run the global SfM pipeline to utilize your algorithm without having to
    implement the rest of the pipeline.

#. What other algorithms can be extended?

    Everything in Theia was built to be modular, so the RANSAC, feature
    extraction, feature matching, pose estimation methods, and more can all be
    easily extended by deriving from the appropriate abstract base classes. Once
    new algorithms are implemented they can be trivially added to Theia and
    automatically be used in the entire pipeline.
