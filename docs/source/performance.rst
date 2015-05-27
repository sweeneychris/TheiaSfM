.. _chapter-performance:

===========
Performance
===========

By utilizing the `Eigen <http://eigen.tuxfamily.org/dox/>`_ and `Ceres Solver <http://www.ceres-solver.org>`_ libraries in addition to custom scalable algorithms, Theia achieves state-of-the-art performance in terms of efficiency and accuracy on large-scale datasets. We measure the performance of Theia in terms of efficiency and accuracy on benchmark datasets for small and large-scale problems to provide meaningful context for the strengths and weaknesses of the library.

.. note:: Results will be posted soon!


Large Scale Benchmarks
======================

We use the `1DSfM Datasets <_http://www.cs.cornell.edu/projects/1dsfm/>`_ as benchmarks for large-scale reconstructions. These datasets provide 2-view matches and epipolar geometry as input, and a reference reconstruction from incremental SfM (computed with `Bundler <http://www.cs.cornell.edu/~snavely/bundler/>`_) for measuring error. The reference reconstruction is not necessarily considered ground truth, but it is a meaningful reference point since incremental SfM algorithms are rather robust and accurate.

We measure the accuracy of camera positions (approximately in meters) and timing results on a 2008 Mac Pro using 1 core (sadly, clang does not allow multi-threading with OpenMP so Ceres cannot utilize multiple cores on this machine). We report N, the number of cameras that could be succesfully reconstructed,in addition to the mean and median camera position errors after robust alignment (via RANSAC) to ground truth camera positions.

.. tabularcolumns:: |l|c||c|c|c|c|

================= ========== === ================ ============== ==========
Dataset           N (input)   N  Median Error (m) Mean Error (m) Timing (s)
================= ========== === ================ ============== ==========
Piccadilly        2152
Union Square      789
Roman Forum       1084
Vienna Cathedral  836
Piazza del Popolo 328
NYC Library       332
Alamo             577
Metropolis        341
Yorkminster       437
Montreal N.D.     450
Tower of London   572
Ellis Island      227
Notre Dame        553
================= ========== === ================ ============== ==========


Small Dataset Benchmarks
========================

To demonstrate the viability of Theia for small scenes, we measure the performance on the `Strecha MVS datasets <http://cvlabwww.epfl.ch/data/multiview/denseMVS.html>`_ datasets. These datasets consist of small scale scenes, though they are extremely high-resolution and rather dense image sampling so it may be considered an "easy" dataset by some. Nonetheless, it is a commonly used dataset to measure reconstruction accuracy.

================= ========== === ================ ============== ==========
Dataset           N (input)   N  Median Error (m) Mean Error (m) Timing (s)
================= ========== === ================ ============== ==========
Castle-19         19
Castle-30         30
Entry-10          10
Fountain-11       11
Herz-Jesu-25      25
Herz-Jesu-8       8
================= ========== === ================ ============== ==========
