.. _chapter-releases:

========
Releases
========

`0.6.0  <https://github.com/sweeneychris/TheiaSfM/archive/v0.6.tar.gz>`_
========================================================================

New Features
------------
* Users can now specify which intrinsics to optimize (using bitmask operators) during bundle adjustment.
* TrackEstimator can now estimate specific tracks (used for incremental SfM).
* Features can be read/written from/to files.
* Matching features can now utilize feature files (so that out-of-core matching can be done).
* Improved SVD efficiency for the 5 point alg.
* 2-view geometric verification takes in the 3D points for BA as input.
* Homographies are now computed during geometric verification
* Homography inlier count is used to choose a better initial pair for Incremental SfM.
* Robust cost functions may now be used for Bundle Adjustment
* Method to estimate a dominant plane from points (by bnuernberger).
* L1 solver now uses the ADMM method. This results in problems that are generally better conditioned and are much faster at scale.
* Reconstructions now can store and compute colors of 3D points.
* Support for reading and writing NVM files (Visual SfM's default format).
* Defaults polynomial solver is now Jenkins Traub method using the `RPolyPlusPlus <https://github.com/sweeneychris/RpolyPlusPlus>`_ implementation.
* Global rotations solvers use a minimum spanning tree to initialize the global orientations.
* Robust cost functions and camera intrinsic optimizations can be used for any type of Bundle Adjustment (not just global BA).
* A new Quadratic Programming module that uses the ADMM method.
* Least Unsquared Deviations position estimation uses the proper QP method for solving.

Bug Fixes
---------
* Several bug fixes for Windows (thanks to Jonas Scheer and others).
* 2-view BA properly holds the camera's extrinsic parameters constant.
* No pointers are used anymore in the cascade hasher. This prevent rarely occuring segfaults.
* Position estimation now fails when view_pairs is empty (thanks to vfragoso).
* Cascade Hasher works properly with zero descriptors (thanks to anfractuosity).
* Root-SIFT handles zero-norm vectors properly.
* Cascade Hasher initializes properly (thanks to anfractuosity).
* Robust rotations solver initializes the sparse matrix from a triplet list. This results in a dramatic speedup for large scale problems.

Misc.
-----
* Incremental SfM will now exit early if no suitable initial pair can be found
* Computing the maximal parallel rigid subgraph only considers views that are in the largest connected component
* RemovesDisconnectedViewPairs now returns the views that were removed
* Only compute HashedImages once per image (for out of core matching)
* Triangulation algorithms return true on success properly
* The point cloud viewer uses the dominant plane detection to set the ground plane for viewing.
* Upgraded EasyEXIF lib to the latest version.
* New CHOLMOD wrapper for sparse cholesky decomposition. This provides a notable performance increase for large scale problems.

`0.5.0  <https://github.com/sweeneychris/TheiaSfM/archive/v0.5.tar.gz>`_
========================================================================

New Features
------------
* Jenkins Traub polynomial root-finding algorithm.
* Cereal library is now used for all I/O.
* Feature matching can now be done in-core or out-of-core.
* Global SfM was completely refactored to be split into RotationEstimator and PositionEstimator classes. This makes implementing new algs straightforward with automatic integration.
* Least unsquared deviations position estimator.
* Linear rotation estimator.
* Extract maximal parallel subgraphs to determine well-constrained positions for estimation.
* Two point algorithm for absolute pose with known vertical direction.
* LMeds (vfragoso).
* Normalized graph cuts (to be used in the future for hiearchical SfM).
* Massively updated flags files for building reconstructions.
* Ability to specify which image pairs to match.

Bug Fixes
---------
* Disable the unit tests for Optimo (thanks to bvanavery).
* Tons of Windows compilation fixes.
* Bundler file I/O fixes (thanks rajvi).
* Fix potential divide by zeros in the RANSAC interface (thanks to klemmster).

Misc.
-----
* Refactoring of the polynomial root-finding algorithms to make the files easier to follow.
* Improved CMake files (thanks to Ceres authors).
* Removed all binary descriptors. This makes the descriptor interfaces much less of a headache.
* Updated VLFeat to the latest version.

`0.4.0 <https://github.com/sweeneychris/TheiaSfM/archive/v0.4.tar.gz>`_
=======================================================================

New Features
------------
* Incremental SfM pipeline.
* New website: `www.theia-sfm.org <http://www.theia-sfm.org>`_.
* Linear method for camera pose registration [JiangICCV]_.
* Better rendering for point clouds.
* Significantly better Cmake scripts for Windows (thanks to bvanevery for testing)
* Mutable priority queue class.
* Bundle adjustment method for cameras only (points held constant).
* Calibrated and Uncalibrated absolute pose estimators.
* Two-view bundle adjustment will now optimize camera intrinsics if they are not known.
* New small and large-scale benchmarking results on the Theia website.

Bug Fixes
---------
* Some Visual Studio bugs and incompatabilities (thanks to Pierre Moulon and Brojeshwar Bhowmick).
* Sample Consensus estimators were incorrectly counting the number of samples needed (found by inspirit).
* Proper normalization the 1dSfM axis of projection.
* OpenGL viewer properly sets zero-values of matrices upon initialization.
* Relative translation optimization (with known rotation) is dramatically improved (thanks to Onur Ozyesil)
* Translations solver uses SPARSE_NORMAL_CHOLESKY when no 3D points are used.

`0.3.0 <https://github.com/sweeneychris/TheiaSfM/archive/v0.3.tar.gz>`_
=======================================================================

New Features
------------
* All cameras are calibrated from EXIF or a median focal length.
* Triangulation is set to use the midpoint method by default.
* All operations on two-view geometry directly operate on the view graph.
* Power method for computing the dominant eigenvector of densor or sparse matrices.
* New program to verify the 1dsfm input against the ground truth model.
* New program to compare two SfM models.
* Nonlinear position estimation uses the nonlinear solver of [WilsonECCV2014]_.
* Removed confusing CameraIntrinsics struct and now all methods use CameraIntrinsicsPrior.
* Calibration files now accept radial distortion and all other camera intrinsics.
* Several new applications to evaluate model and matching quality.
* Robust reconstruction alignment (using RANSAC) to align reconstruction with potential outliers.
* Ability to normalize reconstructions to approximately center and scale nicely for viewing.

Bug Fixes
---------
* 1dSfM dataset input was previously mal-formed.
* GFlags now links pthreads properly.
* Two-view bundle adjustment will no longer use poorly triangulated points for optimization.
* Installation to user-specified folder is done properly.
* Viewing angle test for triangulation.
* Properly estimating relative pose of partially calibrated image matches.

`0.2.0 <https://github.com/sweeneychris/TheiaSfM/archive/v0.2.tar.gz>`_
=======================================================================

New Features
------------

* L1 Solver
* Robust Rotation Solver of [ChatterjeeICCV13]_
* Gflags can now have any namespace
* Reconstructions viewer is now improved
* Initializing rotations from a view graph now use the maximum spanning tree
  instead of a random spanning tree
* Additional run-time options added for building reconstructions

  * ``only_calibrated_views`` will only use calibrated views (from EXIF or
    elsewhere) for building a reconstruction.
  * ``reconstruct_largest_connected_component`` will only build the largest
    connected component of the model instead of building as many models as
    possible.

* 1dSfM datasets [WilsonECCV2014]_ now can be input properly (no quality
  guarantees on the reconstructions though)
* PLY files can be written from a Reconstruction (3D points are all black at
  this point)

Bug Fixes
---------

* Bug fix: removing disconnected view pairs
* Bug fix: 1dSfM filtering of [WilsonECCV2014]_ uses a gaussian distribution to
  randomly sample axis of projections.
* Lowes ratio is fixed.
* Proper hash function for std::pairs (inspiration from Boost)
* Fix BRISK compiler warning for GCC 4.9.1
* Reconstruction viewer bugs and controls are improved
* Better memory management for descriptor extraction and matching

0.1.0
=====

Initial release.
