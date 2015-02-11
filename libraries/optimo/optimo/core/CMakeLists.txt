# Tests:
IF (${OPTIMO_BUILD_TESTING})
  GTEST(objects_tests)
  GTEST(numerical_gradient_tests)
  GTEST(numerical_hessian_tests)
ENDIF (${OPTIMO_BUILD_TESTING})
