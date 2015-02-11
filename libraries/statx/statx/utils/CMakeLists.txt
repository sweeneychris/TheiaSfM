IF (${STATX_BUILD_TESTING})
  GTEST(common_funcs_tests statx)
  GTEST(ecdf_tests statx)
ENDIF (${STATX_BUILD_TESTING})