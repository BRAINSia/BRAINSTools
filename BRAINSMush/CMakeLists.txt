project(BRAINSMush)
#-----------------------------------------------------------------------------
# Output directories.
#

set(ALL_PROGS_LIST
  BRAINSMush
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

if(0)
if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
endif()
