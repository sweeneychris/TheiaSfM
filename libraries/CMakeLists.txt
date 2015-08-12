if (${BUILD_TESTING})
  add_subdirectory(gtest)
endif (${BUILD_TESTING})

# Add Agast keypoint detector
add_subdirectory(agast)

# CImg for image input/output
add_subdirectory(cimg)

# EasyExif for reading image metadata
add_subdirectory(easyexif)

# Add Optimo
add_subdirectory(optimo)

# Add Statx
add_subdirectory(statx)

# STLPlus for filepath tools
add_subdirectory(stlplus3)

# Add VLFeat
add_subdirectory(vlfeat)
