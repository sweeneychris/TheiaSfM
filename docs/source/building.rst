.. _chapter-building:

======================
Building Theia Library
======================

Theia source code and documentation are hosted on `Github
<https://github.com/sweeneychris/TheiaSfM>`_ where you can always grab the latest version

.. _section-dependencies:

Dependencies
------------

Theia relies on a number of open source libraries. Luckily, most of the will be included in Ceres

1. C++11 is needed for certain functionality and added models to the stdlib. C++0x will probably work in most cases, but is not guaranteed. As such, you need a compiler that supports C++11 appropriately.

2. `CMake <http://www.cmake.org>`_ is a cross platform build system. Theia needs a relatively recent version of CMake (version 2.8.0 or better).


3. `eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ is used extensively for doing nearly all the matrix and linear algebra operations.

4. `OpenImageIO <http://www.openimageio.org/>`_ is used to read and write image files.

5. `Ceres Solver <https://code.google.com/p/ceres-solver/>`_ is a library for solving non-linear least squares problems. In particular, Theia uses it for Bundle Adjustment.

**NOTE**: Theia also depends on the following libraries, but they are included in the installation of Ceres so it is likely that you do not need to reinstall them.


6. `google-glog <http://code.google.com/p/google-glog>`_ is used for error checking and logging. Ceres needs glog version 0.3.1 or later. Version 0.3 (which ships with Fedora 16) has a namespace bug which prevents Ceres from building.


7. `gflags <http://code.google.com/p/gflags>`_ is a library for processing command line flags. It is used by some of the examples and tests.


Make sure all of these libraries are installed properly before proceeding. Improperly installing any of these libraries can cause Theia to not build.

.. _section-building:

Building
--------

Building should be equivalent on all platforms, thanks to CMake. To install Theia, simply run the following commands after you have installed the :ref:`section-dependencies`.

First, navigate to the source directory of the Theia library. Then execute the following commands:

.. code-block:: bash

 mkdir theia-build
 cd theia-build
 cmake ..
 make -j4
 make test

If all tests pass, then you are ready to install. Theia can be install using the make install command

.. code-block:: bash

 make install

You can also try running the unit tests individually. The executables should be located in the bin directory of the theia-build folder.


.. _section-customizing:

Customizing the build
---------------------

It is possible to customize the build process by passing appropriate flags to
``CMake``. Use these flags only if you really know what you are doing.


#. ``-DBUILD_TESTING=OFF``: Use this flag to enable or disable building the unit tests. By default, this option is enabled.

#. ``-DBUILD_DOCUMENTATION=ON``: Turn this flag to ``ON`` to build the documentation with Theia. This option is disabled by default.
