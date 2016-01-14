
##- project(ImageCalculator)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(ImageCalculator_ITK
  ITKImageCompare
  ITKImageFilterBase
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKSmoothing
  ITKTestKernel
)
add_executable(ImageCalculator
  ImageCalculator.cxx
  ImageCalculatorProcess2D.cxx
  ImageCalculatorProcess3D.cxx
  ImageCalculatorUtils.cxx)
target_link_libraries(ImageCalculator ${ImageCalculator_ITK_LIBRARIES} )

add_executable(ImageGenerate ImageGenerate.cxx )
target_link_libraries(ImageGenerate ${ImageCalculator_ITK_LIBRARIES})

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
    add_subdirectory(TestSuite)
endif()

