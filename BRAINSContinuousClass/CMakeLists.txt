#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(BRAINSContinuousClass_ITK
  ITKCommon
  ITKIOImageBase
  ITKImageFilterBase
  ITKImageIntensity
  ITKImageStatistics
  ITKStatistics
  ITKTestKernel
  ITKThresholding
)

add_subdirectory(liblinear-1.8)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/liblinear-1.8)

set(ALL_PROGS_LIST
    BRAINSContinuousClass
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES liblinear ${BRAINSContinuousClass_ITK_LIBRARIES} BRAINSCommonLib)
endforeach()

### Add the Testing Subdirectory.
if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()

