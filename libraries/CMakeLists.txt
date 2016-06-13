if (${BUILD_TESTING})
  add_subdirectory(gtest)
endif (${BUILD_TESTING})

# AKAZE feature extractor.
add_subdirectory(akaze)

# Cereal for portable IO.
add_subdirectory(cereal)

# CImg for image input/output.
add_subdirectory(cimg)

# EasyExif for reading image metadata.
add_subdirectory(easyexif)

# Flann for fast approximate nearest neighbor searches.
add_subdirectory(flann)

# Add Optimo.
add_subdirectory(optimo)

# Add Statx.
add_subdirectory(statx)

# STLPlus for filepath tools.
add_subdirectory(stlplus3)

# Add VLFeat.
add_subdirectory(vlfeat)

# Add VisualSfM files.
add_subdirectory(visual_sfm)