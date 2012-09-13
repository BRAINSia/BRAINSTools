project(BRAINSFit)
#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSFitCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSFitCommonLib.h
  )

set(ALL_PROGS_LIST
  BRAINSFit
  BRAINSFitEZ
  )

set(BRAINSFitLibraries BRAINSCommonLib ${ANTS_LIBS})
DebugImageViewerLibAdditions(BRAINSFitLibraries)
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES ${BRAINSFitLibraries})
endforeach()

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
