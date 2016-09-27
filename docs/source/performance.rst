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
<http://theia-sfm.org/build_reconstruction_flags_strecha.txt>`_. Each
reconstruction generated from this file was then compared to the ground truth
reconstruction (which can be generated from the Strecha dataset using the
"create_reconstruction_from_strecha_dataset.cc" program). The camera intrinsic
calibration files can be easily generated from the information provided by
datasets.

.. csv-table:: Strecha Dataset Performance
    :header: Dataset, N (input), N, Median Error (mm), Mean Error (mm), Timing (s)
    :stub-columns: 1

    Castle-19, 19, 19, 14.7, 25.3, 1.54
    Castle-30, 30, 30, 18.5, 21.7, 3.77
    Entry-10, 10, 10, 4.8, 6.0, 1.15
    Fountain-11, 11, 11, 2.0, 2.4, 1.76
    Herz-Jesu-25, 25, 25, 5.1, 5.1, 2.49
    Herz-Jesu-8, 8, 8, 1.9, 3.1, 0.59

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
<http://theia-sfm.org/build_1dsfm_reconstruction_flags.txt>`_. Each
reconstruction generated was then compared to the ground truth reconstruction
provided in the 1dSfM dataset (these are provided as Bundler files, but can be
converted to Theia reconstructions with the "convert_bundle_file.cc"
program).

.. csv-table:: 1DSfM Dataset Position Error
    :header: Dataset, N (input), N, Median Error (m), Mean Error (m)
    :stub-columns: 1

    Alamo, 577, 558, 0.37, 1.62
    Ellis Island, 227, 220, 4.74, 18.38
    Madrid Metropolis, 341, 321, 0.95, 4.10
    Montreal N.D., 450, 448, 0.41, 0.81
    Notre Dame, 553, 540, 0.20, 0.52
    NYC Library, 332, 321, 0.85, 4.91
    Piazza del Popolo, 328, 326, 1.03, 3.91
    Piccadilly, 2152, 2055, 0.72, 2.67
    Roman Forum, 1084, 1045, 2.19, 9.33
    Tower of London, 572, 456, 1.36, 17.38
    Union Square, 789, 720, 4.9, 10.51
    Vienna Cathedral, 836, 797, 2.55, 13.79
    Yorkminster, 437, 414, 1.37, 4.28
    Trafalgar, 5288, 4716, 5.47, 8.39
    Gendarmenmarkt, 733, 657, 10.09, 35.24

.. csv-table:: 1DSfM Dataset Timings (seconds)
    :header: Dataset, N (input), Rotation, Position, BA, Total
    :stub-columns: 1

    Alamo, 577, 3.31, 44.74, 413.05, 497.11
    Ellis Island, 227, 0.51, 4.97, 13.63, 28.34
    Madrid Metropolis, 341, 1.24, 5.75, 33.87, 47.15
    Montreal N.D., 450, 1.42, 23.19, 107.00, 163.82
    Notre Dame, 553, 4.91, 43.37, 196.22, 330.71
    NYC Library, 332, 0.45, 4.35, 46.78, 61.60
    Piazza del Popolo, 328, 0.47, 8.37, 46.30, 61.31
    Piccadilly, 2152, 49.56, 129.21, 72.26, 330.33
    Roman Forum, 1084, 2.03, 23.49, 183.48, 244.41
    Tower of London, 572, 0.47, 8.03, 129.65, 154.45
    Union Square, 789, 1.06, 6.26, 26.82, 47.56
    Vienna Cathedral, 836, 5.46, 41.06, 110.89, 243.83
    Yorkminster, 437, 0.55, 10.17, 59.41, 92.39
    Trafalgar, 5288, 156.331, 387.29, 142.10, 880.74
    Gendarmenmarkt, 733, 1.88, 13.89, 43.32, 72.04
