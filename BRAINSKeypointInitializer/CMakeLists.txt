##- project(BRAINSKeypointInitializer)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSKeypointInitializer_ITK
  ITKCommon
  ITKEigen3
  ITKImageFunction
  ITKImageGrid
  ITKIOImageBase
  ITKIOTransformBase
  ITKTransform
  ITKTransformFactory
)

add_library(BRAINSKeypointLib STATIC
  BRAINSSIFT3DFeatureDetector.cxx
  BRAINSKeypointMatcher.cxx
)
target_include_directories(BRAINSKeypointLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(BRAINSKeypointLib PUBLIC ${BRAINSKeypointInitializer_ITK_LIBRARIES})

StandardBRAINSBuildMacro(NAME BRAINSKeypointInitializer
  TARGET_LIBRARIES BRAINSKeypointLib BRAINSCommonLib
)

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()
