
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

install(TARGETS ImageCalculator ImageGenerate
        RUNTIME DESTINATION ${BRAINSTools_CLI_INSTALL_RUNTIME_DESTINATION} COMPONENT RuntimeLibraries
        LIBRARY DESTINATION ${BRAINSTools_CLI_INSTALL_LIBRARY_DESTINATION} COMPONENT RuntimeLibraries
        ARCHIVE DESTINATION ${BRAINSTools_CLI_INSTALL_ARCHIVE_DESTINATION} COMPONENT Development
        )

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
    add_subdirectory(TestSuite)
endif()

