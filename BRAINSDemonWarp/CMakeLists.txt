project(BRAINSDemonWarp)
#-----------------------------------------------------------------------------
# Output directories.
#

###
option(BRAINSDEMONWARP_DEBUG_OUTPUT "Enable writing out intermediate results in BRAINSDemonWarp" OFF)
configure_file(${BRAINSDemonWarp_SOURCE_DIR}/CTestCustom.ctest ${BRAINSDemonWarp_BINARY_DIR}/CTestCustom.ctest COPYONLY)
configure_file(${BRAINSDemonWarp_SOURCE_DIR}/debugImage.h.in ${BRAINSDemonWarp_BINARY_DIR}/debugImage.h)

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
  endif(NOT VTK_FOUND)
  set_source_files_properties(${BRAINSDemonWarpTemplatesLIB_SRCS} BRAINSDemonWarp.cxx VBRAINSDemonWarp.cxx
    PROPERTIES COMPILE_DEFINITIONS USE_DebugImageViewer=1)
endif()

add_library(BRAINSDemonWarpTemplatesLIB STATIC ${BRAINSDemonWarpTemplatesLIB_SRCS})
if(USE_DebugImageViewer)
  set(VTKLIB ${VTK_LIBRARIES})
endif()
target_link_libraries(BRAINSDemonWarpTemplatesLIB BRAINSCommonLib ${VTK_LIBRARIES})

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
endif(BUILD_TESTING)
