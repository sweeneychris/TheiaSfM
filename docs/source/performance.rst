.. _chapter-performance:

===========
Performance
===========

This page was last updated March 30th, 2016.

By utilizing the `Eigen <http://eigen.tuxfamily.org/dox/>`_ and `Ceres Solver
<http://www.ceres-solver.org>`_ libraries in addition to custom scalable
algorithms, Theia achieves state-of-the-art performance in terms of efficiency
and accuracy on large-scale datasets. We measure the performance of Theia in
terms of efficiency and accuracy on benchmark datasets for small and large-scale
problems to provide meaningful context for the strengths and weaknesses of the
library.


Small Dataset Benchmarks
========================

To demonstrate the viability of Theia for small scenes, we measure the
performance on the `Strecha MVS datasets
<http://cvlabwww.epfl.ch/data/multiview/denseMVS.html>`_ datasets. These
datasets consist of small scale scenes, though they are extremely
high-resolution and rather dense image sampling so it may be considered an
"easy" dataset by some. Nonetheless, it is a commonly used dataset to measure
reconstruction accuracy.

All reconstructions were generated using the parameters found in `this
configuration file
<http://homes.cs.washington.edu/~csweeney/build_reconstruction_flags_strecha.txt>`_. Each
reconstruction generated from this file was then compared to the ground truth
reconstruction (which can be generated from the Strecha dataset using the
"create_reconstruction_from_strecha_dataset.cc" program). The camera intrinsic
calibration files can be easily generated from the information provided by
datasets.

.. csv-table:: Strecha Dataset Performance
    :header: Dataset, N (input), N, Median Error (mm), Mean Error (mm), Timing (s)
    :stub-columns: 1

    Castle-19, 19, 19, 28.2, 36.1, 2.91
    Castle-30, 30, 30, 21.0, 33.5, 4.53
    Entry-10, 10, 10, 4.2, 6.5, 1.76
    Fountain-11, 11, 11, 1.6, 2.6, 2.99
    Herz-Jesu-25, 25, 25, 4.3, 5.4, 4.53
    Herz-Jesu-8, 8, 8, 1.9, 3.6, 1.31

Large Scale Benchmarks
======================

We use the `1DSfM Datasets <http://www.cs.cornell.edu/projects/1dsfm/>`_ as
benchmarks for large-scale reconstructions. These datasets provide 2-view
matches and epipolar geometry as input, and a reference reconstruction from
incremental SfM (computed with `Bundler
<http://www.cs.cornell.edu/~snavely/bundler/>`_) for measuring error. The
reference reconstruction is not necessarily considered ground truth, but it is a
meaningful reference point since incremental SfM algorithms are rather robust
and accurate.

We measure the accuracy of camera positions (approximately in meters) and timing
results. We report N, the number of cameras that could be succesfully
reconstructed, in addition to the mean and median camera position errors after
robust alignment (via RANSAC) to ground truth camera positions.

All reconstructions were generated using the "build_1dsfm_reconstruction.cc"
program and the parameters found in `this config file
<http://homes.cs.washington.edu/~csweeney/build_1dsfm_reconstruction_flags.txt>`_. Each
reconstruction generated was then compared to the ground truth reconstruction
provided in the 1dSfM dataset (these are provided as Bundler files, but can be
converted to Theia reconstructions with the "convert_bundle_file.cc"
program).

.. csv-table:: 1DSfM Dataset Position Error
    :header: Dataset, N (input), N, Median Error (m), Mean Error (m)
    :stub-columns: 1

    Alamo, 577, 561, 0.37, 1.62
    Ellis Island, 227, 220, 1.89, 5.78
    Madrid Metropolis, 341, 321, 1.06, 2.92
    Montreal N.D., 450, 447, 0.42, 0.73
    Notre Dame, 553, 547, 0.22, 0.64
    NYC Library, 332, 321, 0.78, 6.59
    Piazza del Popolo, 328, 328, 1.02, 4.59
    Piccadilly, 2152, 2055, 0.72, 2.67
    Roman Forum, 1084, 1052, 1.04, 7.75
    Tower of London, 572, 448, 1.34, 14.12
    Union Square, 789, 720, 3.18, 8.42
    Vienna Cathedral, 836, 804, 1.99, 9.52
    Yorkminster, 437, 418, 1.39, 4.90
    Trafalgar, 5288, 4788, 4.49, 8.45
    Gendarmenmarkt, 733, 661, 10.69, 32.65

.. csv-table:: 1DSfM Dataset Timings (seconds)
    :header: Dataset, N (input), Rotation, Position, BA, Total
    :stub-columns: 1

    Alamo, 577, 3.35, 45.22, 767.53, 874.41
    Ellis Island, 227, 0.50, 4.59, 70.50, 94.34
    Madrid Metropolis, 341, 1.27, 5.86, 69.62, 95.09
    Montreal N.D., 450, 1.41, 23.87, 128.69, 207.25
    Notre Dame, 553, 4.48, 44.38, 177.73, 372.56
    NYC Library, 332, 0.46, 5.31, 143.60, 194.58
    Piazza del Popolo, 328, 0.48, 8.67, 63.88, 89.19
    Piccadilly, 2152, 49.90, 132.27, 1064.75, 1427.17
    Roman Forum, 1084, 1.96, 23.78, 1196.35, 1302.7
    Tower of London, 572, 0.48, 5.51, 134.72, 201.56
    Union Square, 789, 1.16, 7.65, 73.52, 131.32
    Vienna Cathedral, 836, 5.59, 41.85, 531.60, 764.49
    Yorkminster, 437, 0.61, 10.06, 94.17, 164.33
    Trafalgar, 5288, 155.72, 399.03, 544.58, 1494.19
    Gendarmenmarkt, 733, 1.18, 14.63, 137.07, 202.261
