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

.. note:: These results were computed using only 1 thread. Sadly, clang does not allow multi-threading with OpenMP so Ceres and other libraries cannot utilize multiple cores on the machine that the experiments were run on. We will post results from the multi-threaded experiments once they are available.


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
<http://cs.ucsb.edu/~cmsweeney/build_strecha_reconstructions.txt>`_. Each
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
could be succesfully reconstructed, in addition to the mean and median camera
position errors after robust alignment (via RANSAC) to ground truth camera
positions.

All reconstructions were generated using the "build_1dsfm_reconstruction.cc"
program and the parameters found in `this config file
<http://cs.ucsb.edu/~cmsweeney/build_1dsfm_reconstructions.txt>`_. Each
reconstruction generated was then compared to the ground truth reconstruction
provided in the 1dSfM dataset (these are provided as Bundler files, but can be
converted to Theia reconstructions with the "convert_bundle_file.cc"
program).

.. note:: Results for the first three datasets will be posted soon!

.. tabularcolumns:: |l|c||c|c|c|c|

================= ========== ==== ================ ============== ==========
Dataset           N (input)   N   Median Error (m) Mean Error (m) Timing (s)
================= ========== ==== ================ ============== ==========
Piccadilly        2152
Union Square      789
Roman Forum       1084
Vienna Cathedral  836        821  4.08             13.7           2008
Piazza del Popolo 328        325  1.27             9.2            212
NYC Library       332        325  1.26             2.26           486
Alamo             577        572  0.39             4.44           1060
Metropolis        341        325  4.13             10.7           397
Yorkminster       437        428  2.79             19.5           695
Montreal N.D.     450        442  0.48             3.9            1624
Tower of London   572        458  1.85             17.2           667
Ellis Island      227        218  2.03             7.8            161
Notre Dame        553        538  1.61             7.51           1012
================= ========== ==== ================ ============== ==========
