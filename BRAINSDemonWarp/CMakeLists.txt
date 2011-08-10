project(BRAINSDemonWarp)
set(LOCAL_PROJECT_NAME BRAINSDemonWarp)

#-----------------------------------------------------------------------------
# Output directories.
#

###
configure_file(${BRAINSDemonWarp_SOURCE_DIR}/CTestCustom.ctest ${BRAINSDemonWarp_BINARY_DIR}/CTestCustom.ctest COPYONLY)

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

add_library(BRAINSDemonWarpTemplatesLIB STATIC ${BRAINSDemonWarpTemplatesLIB_SRCS})
target_link_libraries(BRAINSDemonWarpTemplatesLIB ${ITK_LIBRARIES})

if(0) # Build against Slicer
  ## Include the Slicer macro for setting up default locations!
  SlicerMacroBuildCLI(
    NAME BRAINSDemonWarp
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib BRAINSDemonWarpTemplatesLIB ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
  )
  SlicerMacroBuildCLI(
    NAME VBRAINSDemonWarp
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib BRAINSDemonWarpTemplatesLIB ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
  )
else()
  SEMMacroBuildCLI(
    NAME BRAINSDemonWarp
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib BRAINSDemonWarpTemplatesLIB ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
    )
  SEMMacroBuildCLI(
    NAME VBRAINSDemonWarp
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib BRAINSDemonWarpTemplatesLIB ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
    )
endif()

### Add the Testing Subdirectory.
if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif(BUILD_TESTING)

if(NOT Slicer3_SOURCE_DIR)
  add_subdirectory(Example)
endif(NOT Slicer3_SOURCE_DIR)

