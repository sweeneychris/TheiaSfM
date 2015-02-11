# Tests:
IF (${STATX_BUILD_TESTING})
  GTEST(gamma_tests statx)
  GTEST(rayleigh_tests statx)
ENDIF (${STATX_BUILD_TESTING})

#See root directory for macros 
ADD_SUBDIRECTORY(evd)


