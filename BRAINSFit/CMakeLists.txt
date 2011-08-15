project(BRAINSFit)
set(LOCAL_PROJECT_NAME BRAINSFit)
#-----------------------------------------------------------------------------
# Output directories.
#
#SETOPTIONALDEBUGIMAGEVIEWER()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSFitCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSFitCommonLib.h
  )

set(ALL_PROGS_LIST
  BRAINSFit
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
