project(BRAINSROIAuto)
set(LOCAL_PROJECT_NAME BRAINSROIAuto)

#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSROIAutoCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSROIAutoCommonLib.h
  )


set(CLP BRAINSROIAuto)
set(${CLP}_SOURCE BRAINSROIAuto.cxx )
set(${CLP}_PRIMARY_SOURCE BRAINSROIAuto.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSROIAuto.xml "${${CLP}_SOURCE}" "${${CLP}_PRIMARY_SOURCE}" "main" "")
target_link_libraries( ${CLP} BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES})

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
