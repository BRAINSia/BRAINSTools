#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  ITKImageFilterBase
  ITKImageIntensity
  ITKImageStatistics
  ITKStatistics
  ITKThresholding
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})


add_subdirectory(liblinear-1.8)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/liblinear-1.8)

set(ALL_PROGS_LIST
    BRAINSContinuousClass
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES liblinear ${ITK_LIBRARIES} BRAINSCommonLib)
endforeach()

### Add the Testing Subdirectory.
if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()

