.. _chapter-releases:

========
Releases
========

HEAD
====

New Features
------------

* L1 Solver
* Robust Rotation Solver of [ChatterjeeICCV13]_
* Gflags can now have any namespace
* Reconstructions viewer is now improved
* Initializing rotations from a view graph now use the maximum spanning tree
  instead of a random spanning tree
* Additional run-time options added for building reconstructions

Bug Fixes
---------

* Bug fix: removing disconnected view pairs
* Bug fix: 1dSfM filtering of [WilsonECCV2014]_ uses a gaussian distribution to
  randomly sample axis of projections.
* Lowes ratio is fixed.
* Proper hash function for std::pairs (inspiration from Boost)
* Fix BRISK compiler warning for GCC 4.9.1

0.1.0
=====

Initial release.
