project(BRAINSMultiModeSegment)
set(LOCAL_PROJECT_NAME BRAINSMultiModeSegment)

#-----------------------------------------------------------------------------
# Output directories.
#
###
if(0) # Build against Slicer
  ## Include the Slicer macro for setting up default locations!
  SlicerMacroBuildCLI(
    NAME BRAINSMultiModeSegment
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
  )
else()
  SEMMacroBuildCLI(
    NAME BRAINSMultiModeSegment
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
    )
endif()

if(BUILD_TESTING)
  #  add_subdirectory(TestSuite)
endif(BUILD_TESTING)
