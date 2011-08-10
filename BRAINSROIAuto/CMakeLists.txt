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
if(0) # Build against Slicer
  ## Include the Slicer macro for setting up default locations!
  SlicerMacroBuildCLI(
    NAME BRAINSROIAuto
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
  )
else()
  SEMMacroBuildCLI(
    NAME BRAINSROIAuto
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
    )
endif()

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
