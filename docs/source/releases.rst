.. _chapter-releases:

========
Releases
========

HEAD
====


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
