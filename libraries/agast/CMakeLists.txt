project (agast)

INCLUDE_DIRECTORIES(./agast)
SET (agast_sources
  agast/agast5_8.cc
  agast/agast5_8_nms.cc
  agast/agast7_12d.cc
  agast/agast7_12d_nms.cc
  agast/agast7_12s.cc
  agast/agast7_12s_nms.cc
  agast/AstDetector.cc
  agast/nonMaximumSuppression.cc
  agast/oast9_16.cc
  agast/oast9_16_nms.cc
  agast/agast5_8.h
  agast/agast7_12d.h
  agast/agast7_12s.h
  agast/AstDetector.h
  agast/cvWrapper.h
  agast/oast9_16.h)

ADD_LIBRARY(agast ${agast_sources})
INSTALL(TARGETS agast
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
