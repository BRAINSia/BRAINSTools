project(BRAINSMush)
#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.ctest
  ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.ctest COPYONLY)

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
