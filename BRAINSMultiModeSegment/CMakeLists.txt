project(BRAINSMultiModeSegment)
set(LOCAL_PROJECT_NAME BRAINSMultiModeSegment)

#-----------------------------------------------------------------------------
# Output directories.
#

set(CLP BRAINSMultiModeSegment)
set(${CLP}_SOURCE BRAINSMultiModeSegment.cxx )
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSMultiModeSegment.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries( ${CLP} BRAINSCommonLib ${OPTIONAL_DEBUG_LINK_LIBRARIES})

if(BUILD_TESTING)
  #  add_subdirectory(TestSuite)
endif(BUILD_TESTING)
