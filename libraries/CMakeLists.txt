IF (${BUILD_TESTING})
  ADD_SUBDIRECTORY(gtest)
ENDIF (${BUILD_TESTING})

# Add Agast keypoint detector
ADD_SUBDIRECTORY(agast)

# CImg for image input/output
ADD_SUBDIRECTORY(cimg)

# EasyExif for reading image metadata
ADD_SUBDIRECTORY(easyexif)

# Add Optimo
ADD_SUBDIRECTORY(optimo)

# Add Statx
ADD_SUBDIRECTORY(statx)

# STLPlus for filepath tools
ADD_SUBDIRECTORY(stlplus3)

# Add VLFeat
ADD_SUBDIRECTORY(vlfeat)
