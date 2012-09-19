##- project(BRAINSDemonWarp)
#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSDemonWarpCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSDemonWarpCommonLib.h
  )


set(BRAINSDemonWarpTemplatesLIB_SRCS
         BRAINSDemonWarpTemplates_Vdouble.cxx
         BRAINSDemonWarpTemplates_Vshort.cxx
         BRAINSDemonWarpTemplates_Vushort.cxx
         BRAINSDemonWarpTemplates_int.cxx
         BRAINSDemonWarpTemplates_uint.cxx
         BRAINSDemonWarpTemplates_Vfloat.cxx
         BRAINSDemonWarpTemplates_Vuchar.cxx
         BRAINSDemonWarpTemplates_double.cxx
         BRAINSDemonWarpTemplates_short.cxx
         BRAINSDemonWarpTemplates_ushort.cxx
         BRAINSDemonWarpTemplates_Vint.cxx
         BRAINSDemonWarpTemplates_Vuint.cxx
         BRAINSDemonWarpTemplates_float.cxx
         BRAINSDemonWarpTemplates_uchar.cxx
         BRAINSDemonWarpTemplates.cxx
)

if(USE_DebugImageViewer)
  if(NOT VTK_FOUND) #See if it has already been found because this is a nested project
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
  endif()
endif()

add_library(BRAINSDemonWarpTemplatesLIB STATIC ${BRAINSDemonWarpTemplatesLIB_SRCS})

set(BRAINSDemonWarpTemplatesLIBLibraries BRAINSCommonLib)
DebugImageViewerLibAdditions(BRAINSDemonWarpTemplatesLIBLibraries)

target_link_libraries(BRAINSDemonWarpTemplatesLIB ${BRAINSDemonWarpTemplatesLIBLibraries})

set(ALL_PROGS_LIST
  BRAINSDemonWarp
  VBRAINSDemonWarp
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSDemonWarpTemplatesLIB )
endforeach()

### Add the Testing Subdirectory.
if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif()
