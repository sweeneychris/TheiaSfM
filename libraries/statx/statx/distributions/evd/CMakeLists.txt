# Tests:
IF (${STATX_BUILD_TESTING})
  GTEST(gpd_tests statx)
  GTEST(gev_tests statx)
  GTEST(gumbel_tests statx)
ENDIF (${STATX_BUILD_TESTING})