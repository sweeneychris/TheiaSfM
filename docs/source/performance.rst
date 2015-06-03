.. _chapter-performance:

===========
Performance
===========

By utilizing the `Eigen <http://eigen.tuxfamily.org/dox/>`_ and `Ceres Solver
<http://www.ceres-solver.org>`_ libraries in addition to custom scalable
algorithms, Theia achieves state-of-the-art performance in terms of efficiency
and accuracy on large-scale datasets. We measure the performance of Theia in
terms of efficiency and accuracy on benchmark datasets for small and large-scale
problems to provide meaningful context for the strengths and weaknesses of the
library.

.. note:: These results were computed using only 1 thread. Sadly, clang does not allow multi-threading with OpenMP so Ceres cannot utilize multiple cores on the machine that the experiments were run on. We will post results from the multi-threaded experiments once they are available.


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
results on a 2008 Mac Pro using 1 core. We report N, the number of cameras that
could be succesfully reconstructed,in addition to the mean and median camera
position errors after robust alignment (via RANSAC) to ground truth camera
positions.

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

To demonstrate the viability of Theia for small scenes, we measure the
performance on the `Strecha MVS datasets
<http://cvlabwww.epfl.ch/data/multiview/denseMVS.html>`_ datasets. These
datasets consist of small scale scenes, though they are extremely
high-resolution and rather dense image sampling so it may be considered an
"easy" dataset by some. Nonetheless, it is a commonly used dataset to measure
reconstruction accuracy.

All reconstructions were generated using the parameters found in `this
configuration file
<http://cs.ucsb.edu/~cmsweeney/build_strecha_reconstructions.txt>`_. The
reconstruction generated from this file was then compared to the ground truth
reconstruction (which can be generated from the Strecha dataset using the
"create_reconstruction_from_strecha_dataset.cc" program). The camera intrinsic
calibration files can be easily generated from the information provided by
datasets.

================= ========== === ================= =============== ==========
Dataset           N (input)   N  Median Error (mm) Mean Error (mm) Timing (s)
================= ========== === ================= =============== ==========
Castle-19         19         19  33.01             38.26           6.51
Castle-30         30         30  29.88             32.47           11.61
Entry-10          10         10  4.91              5.76            3.84
Fountain-11       11         11  2.11              2.43            6.04
Herz-Jesu-25      25         25  5.21              5.30            10.32
Herz-Jesu-8       8          8   3.32              3.58            2.06
================= ========== === ================= =============== ==========
