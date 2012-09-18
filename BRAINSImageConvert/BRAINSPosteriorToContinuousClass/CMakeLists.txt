
StandardBRAINSBuildMacro(NAME BRAINSPosteriorToContinuousClass
TARGET_LIBRARIES ${ITK_LIBRARIES})
#CLI_XML_FILE BRAINSPosteriorToContinuousClass.xml)

if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif()
