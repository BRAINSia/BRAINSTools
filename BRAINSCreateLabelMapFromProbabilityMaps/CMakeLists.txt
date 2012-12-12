StandardBRAINSBuildMacro( NAME BRAINSCreateLabelMapFromProbabilityMaps
  TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES})

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()
