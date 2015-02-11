Copyright (C) 2014 Victor Fragoso <vfragoso@cs.ucsb.edu>
--------------------------------------------------------------------------

STATX
-----

This library implements several MLE estimators for various statistical 
distributions as well as useful statistical functions.

Installation
------------

To build and install this library CMake (http://www.cmake.org) is required.
The following steps can be used to build the library:

1. Invoke Cmake
   $ cmake -DCMAKE_INSTALL_PREFIX=<DESTINATION> .

The <DESTINATION> is the installation path, e.g., /opt/local/.

2. Compile
   $ make

3. Make install
   $ make install

Make sure that you have the right permissions to write in the DESTINATION
directory for the installation.

Dependencies
------------

1. Eigen: http://eigen.tuxfamily.org/

2. Google Flags: https://code.google.com/p/gflags/

3. Google Glog: https://code.google.com/p/google-glog/

4. Optimo: https://bitbucket.org/vfragoso/optimo


Optional dependencies
---------------------

1. ceres-solver: https://code.google.com/p/ceres-solver/

Unit Tests
----------

1. Invoke CMake:
$ cmake -DBUILD_TESTING=TRUE .

2. Invoke Makefile:
$ make

3. Run tests
$ ./bin/<component>_tests

Contact: vfragoso@cs.ucsb.edu

RELEASE NOTES:

- The parameter estimation for the generalized Pareto distribution is not 
yet implemented.