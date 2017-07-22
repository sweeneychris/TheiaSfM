## [Unreleased]
### Added
- Added virtual destructors to the `SymEigsSolver` and `UpperHessenbergQR` classes
  to fix compiler warnings, by [Julian Kent](https://github.com/jkflying)
- Added a `NUMERICAL_ISSUE` entry to the `COMPUTATION_INFO` enumeration to indicate
  the status of Cholesky decomposition
- Added the `info()` member function to `DenseCholesky` and `SparseCholesky` to
  report the status of the decomposition

### Changed
- Documentation updates
- Update project URL to be [https://spectralib.org](https://spectralib.org).


## [0.5.0] - 2017-02-05
### Added
- Added the generalized eigen solver `SymGEigsSolver` in the regular inverse mode
- Added the wrapper class `SparseRegularInverse` that can be used with
  `SymGEigsSolver` in the regular inverse mode
- Added test code for generalized eigen solver in the regular inverse mode

### Changed
- Improved the numerical precision and stability of some internal linear
  algebra classes, including `TridiagEigen`, `UpperHessenbergEigen`, and
  `DoubleShiftQR`
- **API change**: The `x_in` argument in matrix operation functions, e.g.
  `perform_op()`, is now labelled to be constant
- Fixed a [bug](https://github.com/yixuan/spectra/issues/15) that
  `GenEigsComplexShiftSolver` gave wrong results when transforming back the
  eigenvalues, discovered by [@jdbancal](https://github.com/jdbancal)
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.7.0
- Documentation improvement


## [0.4.0] - 2016-11-14
### Added
- Added an `Uplo` template parameter to the `DenseSymShiftSolve` class
- Added the generalized eigen solver `SymGEigsSolver` in the Cholesky mode
- Added the wrapper classes `DenseCholesky` and `SparseCholesky` that can be
  used with `SymGEigsSolver` in the Cholesky mode
- Added test code for generalized eigen solver in the Cholesky mode

### Changed
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.5.7
- Improved documentation
- Updated Travis CI script
- Allowing basic math functions such as `abs()` and `sqrt()` to be overloaded
  (avoid using `std::abs` and `std::sqrt` directly), thanks to
  [@jdbancal](https://github.com/jdbancal). This makes it possible to use
  user-defined float number types with Spectra
- Replaced other `std` functions by their Eigen counterparts, for example using
  `Eigen::NumTraits<Scalar>::epsilon()` to substitute
  `std::numeric_limits<Scalar>::epsilon()`
- Improved the numerical stability of several operations, e.g. the function
  `hypot(x, y)` is used to compute `sqrt(x^2 + y^2)`
- More careful use of "approximate zero" constants
- Fixed an out-of-bound [bug](https://github.com/yixuan/spectra/issues/14)
  detected by [@jdbancal](https://github.com/jdbancal)


## [0.3.0] - 2016-07-03
### Added
- Added the wrapper classes `SparseSymMatProd` and `SparseSymShiftSolve`
  for sparse symmetric matrices
- Added the wrapper class `SparseGenRealShiftSolve` for general sparse matrices
- Added tests for sparse matrices
- Using Travis CI for automatic unit test

### Changed
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.5.6
- **API change**: Each eigen solver was moved to its own header file.
  For example to use `SymEigsShiftSolver` one needs to include
  `<SymEigsShiftSolver.h>`
- Header files for internal use were relocated


## [0.2.0] - 2016-02-28
### Added
- Benchmark script now outputs number of matrix operations
- Added this change log
- Added a simple built-in random number generator, so that the algorithm
  was made to be deterministic
- Added the wrapper class `DenseSymMatProd` for symmetric matrices

### Changed
- Improved Arnoldi factorization
  - Iteratively corrects orthogonality
  - Creates new residual vector when invariant subspace is found
  - Stability for matrices with repeated eigenvalues is greatly improved
- Adjusted deflation tolerance in double shift QR
- Updated result analyzer
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.3.4
- Updated copyright information
- **API change**: Default operator of `SymEigsSolver` was changed from
  `DenseGenMatProd` to `DenseSymMatProd`


## [0.1.0] - 2015-12-19
### Added
- Initial release of Spectra
