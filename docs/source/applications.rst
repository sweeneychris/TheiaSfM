.. _chapter-applications:

============
Applications
============

There are several applications that come with Theia right out of the box. These
applications are useful on their own, and also help provide context for how
Theia can be used for your own applications.

Features
========

Extract Features
----------------

Extract any type of feature that is implemented in Theia (SIFT, BRIEF, BRISK,
FREAK) and display the feature matches by outputting the images and
features. The feature type and number of threads can be set at runtime.

.. code-block:: bash

  ./bin/extract_features --input_imgs=/path/to/images/*.jpg --img_output_dir=/path/to/output --num_threads=4 --descriptor=SIFT

Match Descriptors
-----------------

Given an input set of images, this program will match descriptors between all
images. Many parameters can be set at runtime, and the matched images are
written to a specified output directory.

.. code-block:: bash

  ./bin/match_features --input_imgs=/path/to/images/*.jpg --img_output_dir=/path/to/output --num_threads=4 --descriptor=SIFT --matcher=brute_force --lowes_ratio=0.8

Reconstructions
===============

Build Reconstruction
--------------------

This application will build a 3D reconstruction from a set of images or a set of
image matches. Detailed documentation for the structure-from-motion pipeline can
be found at :ref:`chapter-sfm`. Many parameters can be set at runtime (too many
to list here), and we provide an example of the possible settings in
applications/build_reconstruction_flags.txt. This flags file may be run by
executing the command:

.. code-block:: bash

  ./bin/build_reconstruction --flagfile=/path/to/build_reconstructions_flags.txt

The reconstruction parameters may need to be tuned a bit

If images are supplied as input, then features are extracted and matched between
images before the reconstruction process begins. It is advised that you save
these matches (by specifying a --output_matches_file=/path/to/output.matches) so
that the reconstruction process may be restarted directly from the two-view
geometry. This allows you to tune the reconstruction parameters without having
to wait for image matching which is typically the slowest part of
structure-from-motion. Alternatively, you could first generate the two view
geometry and save the information using the program below.

Compute Two View Geometry
-------------------------

Computes the two view matches and geometry between image pairs. This program
follows many of the same parameters as the Build Reconstructions program, but is
useful for generating two view geometries prior to building a
reconstruction. Feature matching is performed between images then geometric
verification is performed to determine which feature matches are inliers. Only
image pairs that have sufficiently many geometrically-verified matches are
considered valid.

Compute Reconstruction Statistics
---------------------------------

Computes reprojection error, etc.

View Reconstruction
-------------------

A very, very basic OpenGL point cloud viewer. Improvements are very welcome!

Create Calibration File From EXIF
---------------------------------

Creates a calibration file from the EXIF information that can be
extracted from an image set.

Convert Bundle File
-------------------

Converts a bundler reconstruction to a Theia :class:`Reconstruction`.

Convert Sift Key File
---------------------

Converts Lowe's SIFT key files to a binary format.
