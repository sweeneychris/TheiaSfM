.. highlight:: c++

.. default-domain:: cpp

.. _documentation-math:

====
Math
====

At the root of computer vision is a heavy amount of math and probability. Theia contains various math functions implemented with a generic interface for ease of use.

.. _section-closed_form_poly:

Closed Form Polynomial Solver
=============================

Many problems in vision rely on solving a polynomial quickly. For small degrees
(n <= 4) this can be done in closed form, making them exceptionally fast. We
have implemented solvers for these cases.

.. function:: int SolveQuadraticReals(double a, double b, double c, double* roots)

.. function:: int SolveQuadratic(double a, double b, double c, std::complex<double>* roots)

  Provides solutions to the equation :math:`a*x^2 + b*x + c = 0`

.. note:: The closed form solutions for cubic and quartic solvers are known to be numerically unstable. We recommend using the generic polynomial solvers (below) instead. This will sacrifice efficiency a small amount for a significant improvement in solution quality.

.. function::  int SolveCubicReals(double a, double b, double c, double d, double* roots)

.. function::  int SolveCubic(double a, double b, double c, double d, std::complex<double>* roots)

   Provides solutions to the equation :math:`a*x^3 + b*x^2 + c*x + d = 0` using `Cardano's <http://en.wikipedia.org/wiki/Cubic_function#Cardano.27s_method>`_ method.


.. function::  int SolveQuarticReals(double a, double b, double c, double d, double e, double* roots)

.. function::  int SolveQuartic(double a, double b, double c, double d, double e, std::complex<double>* roots)

  Provides solutions to the equation :math:`a*x^4 + b*x^3 + c*x^2 + d*x + e = 0` using `Ferrari's method <http://en.wikipedia.org/wiki/Quartic_function#Ferrari.27s_solution>`_ to reduce to problem to a depressed cubic.


.. _section-generic_poly:

Generic Polynomial Solver
=========================

For polynomials of degree > 4 there are no closed-form solutions, making the
problem of finding roots much more difficult. However, we have implemented
several functions that will solve for polynomial roots. For all polynomials we
require that the largest degree appears first and the smallest degree appears
last in the input VectorXd such that:

.. math:: \sum_{i=0}^N p(i) x^{N-i}

Where :math:`p(i)` is the input VectorXd.

.. function:: bool FindPolynomialRoots(const Eigen::VectorXd& polynomial, Eigen::VectorXd* real, Eigen::VectorXd* imaginary)

  This function finds the roots of the input polynomial using one of the methods
  below. All methods in Theia that require finding polynomial roots use this
  method. This is so that we can easily change the default root-finding method
  of choice (i.e. Companion Matrix to Jenkins-Traub, etc.) by modifying this
  function once instead of modify every instance where we want to find
  polynomial roots. This allows us to easily swap in new polynomial root-solvers
  (that may be more efficient or numerically stable) as they are implemented.

.. function:: bool FindPolynomialRootsJenkinsTraub(const Eigen::VectorXd& polynomial, Eigen::VectorXd* real, Eigen::VectorXd* imaginary)

  The `Jenkins Traub algorithm <https://en.wikipedia.org/wiki/Companion_matrix>`_
  is a three-stage algorithm for finding roots of polynomials with real
  coefficients as outlined in [JenkinsTraub]_. Please note that
  this variant is different than the complex-coefficient version, and is
  estimated to be up to 4 times faster.

  The algorithm works by computing shifts in so-called "K-polynomials" that
  deflate the polynomial to reveal the roots. Once a root is found (or in the
  real-polynomial case, a pair of roots) then it is divided from the polynomial
  and the process is repeated. This method is consider to be "pratically a
  standard in black-box polynomial root finder" (Numerical Recipes 2007) and is
  based on the `Rpoly++ <http://github.com/sweeneychris/RpolyPlusPlus>`_ implementation.

.. function:: bool FindPolynomialRootsCompanionMatrix(const Eigen::VectorXd& polynomial, Eigen::VectorXd* real, Eigen::VectorXd* imaginary)

  Roots are computed using the `Companion Matrix <https://en.wikipedia.org/wiki/Companion_matrix>`_ with balancing to help improve
  the condition of the matrix system we solve. This is a reliable, stable method
  for computing roots but is most often the slowest method.

.. function:: double FindRootIterativeLaguere(const Eigen::VectorXd& polynomial, const double x0, const double epsilon, const int max_iter)

  Finds a single polynomials root iteratively based on the starting position :math:`x_0` and
  guaranteed precision of epsilon using `Laguerre's Method <https://en.wikipedia.org/wiki/Laguerre%27s_method>`_.

.. function:: double FindRootIterativeNewton(const Eigen::VectorXd& polynomial, const double x0, const double epsilon, const int max_iter)

  Finds a single polynomials root iteratively based on the starting position :math:`x_0` and
  guaranteed precision of epsilon using `Newton's Method <https://en.wikipedia.org/wiki/Newton%27s_method>`_.

.. _section-matrix_methods:

Matrix Methods
==============

Theia implements many useful linear algebra methods including optimizations, factorizations, and utility methods.

.. class:: L1Solver

  We implement a robust :math:`L_1` solver that minimizes :math:`||Ax - b||_1`
  under :math:`L_1` norm. This problem may be cast as a simple and efficient
  linear program and solved with interior point methods. The interface is fairly
  generic and may be used with sparse or dense matrices. An initial guess is
  needed for :math:`x` to perform the minimization.

.. member:: double L1Solver::Options::max_num_iterations

  DEFAULT: ``100``

  The maximum number of iterations to perform before stopping.

.. code:: c++

  Eigen::MatrixXd A;
  Eigen::VectorXd b, x;
  // Fill A and b with known values.

  L1Solver::Options option;
  L1Solver<Eigen::MatrixXd> l1_solver(options, A);
  l1_solver.Solve(b, &x);
  // x now contains the solution that minimizes ||Ax - b|| under L1 norm.


.. class:: DominantEigensolver

  This class finds the dominant eigenvalue/eigenvector pair of a given matrix
  using `power iterations <https://en.wikipedia.org/wiki/Power_iteration>`_. We
  use a generic interface that utilizes the :class:`LinearOperator` class so
  that the user may determine who the dominant eigenvalues are computed. For
  instance, by passing the :class:`SparseInveseLULinearOperator` the
  :class:`DominantEigensolver` performs inverse power iterations and thus the
  smallest eigenvalue/eigenvector pair may be computed. This is useful for
  recovering a null space vector.

.. function:: void GaussJordan(Eigen::MatrixBase<Derived>* input, int max_rows = 99999)

  Perform traditional Gauss-Jordan elimination on an Eigen3 matrix. If
  ``max_rows`` is specified, it will on perform Gauss-Jordan on the first
  ``max_rows`` number of rows. This is useful for problems where your system is
  extremely overdetermined and you do not need all rows to be solved.


.. _section-sprt:

Sequential Probability Ratio Test
=================================

Modified version of Wald's `SPRT <http://en.wikipedia.org/wiki/Sequential_probability_ratio_test>`_ as [Matas]_ et. al. implement it in "Randomized
RANSAC with Sequential Probability Ratio Test"

.. function:: double CalculateSPRTDecisionThreshold(double sigma, double epsilon, double time_compute_model_ratio = 200.0, int num_models_verified = 1)

 ``sigma``: Probability of rejecting a good model (Bernoulli parameter).

 ``epsilon``: Inlier ratio.

 ``time_compute_model_ratio``: Computing the model parameters from a sample takes the same time as verification of time_compute_model_ratio data points. Matas et. al. use 200.

 ``num_model_verified``: Number of models that are verified per sample.

 ``Returns``:  The SPRT decision threshold based on the input parameters.


.. function:: bool SequentialProbabilityRatioTest(const std::vector<double>& residuals, double error_thresh, double sigma, double epsilon, double decision_threshold, int* num_tested_points, double* observed_inlier_ratio)

 Modified version of Wald's SPRT as [Matas]_ et. al. implement it in "Randomized
 RANSAC with Sequential Probability Ratio Test". See the paper for more
 details.

 ``residuals``: Error residuals to use for SPRT analysis.

 ``error_thresh``: Error threshold for determining when Datum fits the model.

 ``sigma``: Probability of rejecting a good model.

 ``epsilon``: Inlier ratio.

 ``decision_threshold``: The decision threshold at which to terminate.

 ``observed_inlier_ratio``: Output parameter of inlier ratio tested.
