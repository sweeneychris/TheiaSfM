Copyright (C) 2014 Victor Fragoso <vfragoso@cs.ucsb.edu>
--------------------------------------------------------------------------

Optimo
------

This header library implements several optimization solvers such as
Newton, Gradient Descent, BFGS, and some primal dual newton solvers for LPs
and QPs. 

Installation
------------

There are two options to install Optimo:

1. Copy the headers preserving the structure directory in optimo/ to some desired 
   installation directory.

2. If you have cmake (easiest way):

   i) Invoke CMake:
      $ cmake -DCMAKE_INSTALL_PREFIX=<DESTINATION> .

   Here <DESTINATION> is the installation path, e.g., /usr/local/include/.

   ii). Install:
      $ make install

   Make sure you have permissions to write to the destination directory.

If you want to remove optimo, then simply delete the <DESTINATION>/optimo directory.

Dependencies
------------

1. Eigen 3: http://eigen.tuxfamily.org/

Optional dependencies:

If interested in using the primal dual LP solvers, then you require the 
following libraries from SuiteSparse:

1. SuiteSparse: https://www.cise.ufl.edu/research/sparse/SuiteSparse/

2. Cholmod: https://www.cise.ufl.edu/research/sparse/cholmod/

3. Umfpack: https://www.cise.ufl.edu/research/sparse/umfpack/

Unit tests
----------

To build this the testing examples this library requires:

a) CMake: http://www.cmake.org

b) Google Flags: https://code.google.com/p/gflags/

c) Google Glog: https://code.google.com/p/google-glog/

Building the Unit Tests

1. Invoke CMake:
$ cmake -DBUILD_TESTING=TRUE .

2. Invoke Makefile:
$ make

3. Run tests
$ ./bin/<component>_tests

Contact: vfragoso@cs.ucsb.edu

RELEASE NOTES: 

- The Infeasible Newton solver needs a more rigorous testing.

