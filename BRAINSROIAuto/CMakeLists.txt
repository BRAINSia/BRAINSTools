project(BRAINSROIAuto)
set(LOCAL_PROJECT_NAME BRAINSROIAuto)

#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSROIAutoCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSROIAutoCommonLib.h
  )

###
set(ALL_PROGS_LIST
  BRAINSROIAuto
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
