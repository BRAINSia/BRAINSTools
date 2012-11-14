
StandardBRAINSBuildMacro(NAME BRAINSPosteriorToContinuousClass
TARGET_LIBRARIES ${ITK_LIBRARIES})
#CLI_XML_FILE BRAINSPosteriorToContinuousClass.xml)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
