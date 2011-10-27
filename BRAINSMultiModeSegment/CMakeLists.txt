project(BRAINSMultiModeSegment)
#-----------------------------------------------------------------------------
# Output directories.
#
###
set(ALL_PROGS_LIST
  BRAINSMultiModeSegment
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

if(BUILD_TESTING)
  #  add_subdirectory(TestSuite)
endif(BUILD_TESTING)
