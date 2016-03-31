.. Theia documentation master file, created by
   sphinx-quickstart on Wed Mar  6 17:39:11 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====
Theia
=====

Theia is a computer vision library developed by `Chris Sweeney <http://cs.ucsb.edu/~cmsweeney>`_ aimed at providing efficient and reliable
algorithms for Structure from Motion (SfM). The goal of this library is to provide researchers with an out of the box tool for multi-view reconstruction that can be easily extended. Many common algorithms for pose, feature detection and description, matching, and reconstruction have been implemented. All contain simple interfaces, limited dependencies, and extensive documentation.

* Download the latest `release <https://github.com/sweeneychris/TheiaSfM>`_ or clone the
  `repo <https://github.com/sweeneychris/Theia>`_ for development.

* Read the :ref:`chapter-api` and the :ref:`chapter-applications` to learn more about Theia

* If you have questions, please email the `Theia mailing list <http://groups.google.com/group/theia-vision-library>`_.

.. _section-Documentation:

=============
Documentation
=============

To use Theia, simply add the following line to your program after you build and
link the library:

``#include <theia/theia.h>``

After compiling and executing ``make install``, the simplest way to add and link
the library is with the following CMake commands:

.. code::

  find_package(Theia REQUIRED)
  include_directories(${THEIA_INCLUDE_DIRS})
  target_link_libraries(my_library ${THEIA_LIBRARIES})

We attempt to provide sufficient documentation but often further documentation
can be found in the source code itself. You will likely find the API
documentation useful as well. Additionally, (nearly) every file is covered by a
unit test that can be viewed as an example use case of the various methods and
classes in Theia. If you have looked at the documentation, the tutorials, the
source code, and the unit tests and still have confusion please email `the Theia
mailing list <http://groups.google.com/group/theia-vision-library>`_

Finally, it should be noted that all the code in Theia is under the namespace
theia, so you will have to reference that namespace in order to use functions
from this library.

Performance
===========

Theia achieves state-of-the-art SfM performance on large-scale
datasets. Efficiency and robustness is a key component of the library. You can
see the latest performance benchmarks for small and large-scale datasets on the
:ref:`chapter-performance` page.

Citation
========

If you use Theia for an academic publication, please cite this
manual. e.g., ::

  @misc{theia-manual,
    author = {Chris Sweeney},
    title = {Theia Multiview Geometry Library: Tutorial \& Reference},
    howpublished = "\url{http://theia-sfm.org}",
  }

When using specific algorithms that are implemented within Theia, we ask that
you please cite the original sources. More information on which files use which
references in the literature can be found in the header files of the relevant
functions.

================
Acknowledgements
================

Theia was originally developed to provide a centralized code base to the `Four
Eyes Lab <http://ilab.cs.ucsb.edu>`_ at UC Santa Barbara, but has since been
expanded to an open-source project for the vision community.

The core of the original library is written by `Chris Sweeney
<http://cs.ucsb.edu/~cmsweeney>`_. Funding for Theia was provided by his
advisors `Tobias Hollerer <http://cs.ucsb.edu/~holl>`_ and `Matthew Turk
<http://cs.ucsb.edu/~mturk>`_ in part by NSF Grant IIS-1219261, ONR Grant
N00014-14-1-0133, and NSF Graduate Research Fellowship Grant DGE-1144085.

.. toctree::
   :maxdepth: 3
   :hidden:

   building
   api
   applications
   faq
   performance
   contributions
   releases
   bibliography
   license
