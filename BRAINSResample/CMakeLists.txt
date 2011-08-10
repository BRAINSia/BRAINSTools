project(BRAINSResample)
set(LOCAL_PROJECT_NAME BRAINSResample)

#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSResampleCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSResampleCommonLib.h
  )

###
configure_file(${BRAINSResample_SOURCE_DIR}/CTestCustom.ctest ${BRAINSResample_BINARY_DIR}/CTestCustom.ctest COPYONLY)

include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBRAINS3BuildMacros.cmake)
###


set(CLP BRAINSResample)
#set(${CLP}_SOURCE BRAINSResample.cxx ${BRAINSResample_SOURCE_DIR}/BRAINSResamplePrimary.cxx)
set(${CLP}_SOURCE BRAINSResample.cxx )
set(${CLP}_PRIMARY_SOURCE ${BRAINSResample_SOURCE_DIR}/BRAINSResamplePrimary.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSResamplePrimary.xml "${${CLP}_SOURCE}" "${${CLP}_PRIMARY_SOURCE}" "BRAINSResamplePrimary" "")

### Add the Testing Subdirectory.
if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif(BUILD_TESTING)
