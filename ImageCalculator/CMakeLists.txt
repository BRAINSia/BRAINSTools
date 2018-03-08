
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

set(MODULE_FOLDER "Module-ImageCalculator")

add_executable(ImageCalculator
  ImageCalculator.cxx
  ImageCalculatorProcess2D.cxx
  ImageCalculatorProcess3D.cxx
  ImageCalculatorUtils.cxx)
target_link_libraries(ImageCalculator ${ImageCalculator_ITK_LIBRARIES} )
set_target_properties(ImageCalculator PROPERTIES FOLDER ${MODULE_FOLDER})

add_executable(ImageGenerate ImageGenerate.cxx )
target_link_libraries(ImageGenerate ${ImageCalculator_ITK_LIBRARIES})
set_target_properties(ImageGenerate PROPERTIES FOLDER ${MODULE_FOLDER})

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
    add_subdirectory(TestSuite)
endif()

