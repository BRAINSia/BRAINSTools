##- project(BRAINSLandmarkInitializer)
#-----------------------------------------------------------------------------
# Output directories.
#
###
set(ALL_PROGS_LIST
  BRAINSLandmarkInitializer
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

#if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
#    add_subdirectory(TestSuite)
#endif()
