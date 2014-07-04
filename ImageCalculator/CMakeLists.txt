
##- project(ImageCalculator)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKCommon
  ITKDeprecated
  ITKIOImageBase
  ITKImageCompare
  ITKImageFilterBase
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKSmoothing
  ITKTestKernel
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

add_executable(ImageCalculator
  ImageCalculator.cxx
  ImageCalculatorProcess2D.cxx
  ImageCalculatorProcess3D.cxx
  ImageCalculatorUtils.cxx)
target_link_libraries(ImageCalculator ${ITK_LIBRARIES} )

add_executable(ImageGenerate ImageGenerate.cxx )
target_link_libraries(ImageGenerate ${ITK_LIBRARIES})

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()

