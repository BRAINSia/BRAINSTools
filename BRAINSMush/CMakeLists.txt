project(BRAINSMush)
set(LOCAL_PROJECT_NAME BRAINSMush)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

enable_testing()
#include(Dart)
include(CTest)
include(CPack)
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBuildMacros.cmake)

###
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

###
CHECKIFSLICER3BUILD()  ## Call the convenience macro

if(NOT ITK_FOUND)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

#-----------------------------------------------------------------------------
# Output directories.
#
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBRAINS3BuildMacros.cmake)
#SETOPTIONALDEBUGIMAGEVIEWER()


# No need to add with nothing to compile. add_subdirectory(BRAINSMush_Common)
include(${BRAINS_CMAKE_HELPER_DIR}/IJMacros.txt)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.ctest ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.ctest COPYONLY)

add_library(BRAINSMushCOMMONLIB BRAINSMushPrimary.cxx BRAINSMushPrimaryCLP.h)
target_link_libraries(BRAINSMushCOMMONLIB BRAINSCommonLib ${ITK_LIBRARIES})

set(CLP BRAINSMush)
set(${CLP}_SOURCE BRAINSMush.cxx )
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSMushPrimary.xml "${${CLP}_SOURCE}" "" "main" "BRAINSMushCOMMONLIB")
target_link_libraries( ${CLP} BRAINSMushCOMMONLIB BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES})

###
add_executable(BRAINSMushTest BRAINSMushTest.cxx)
target_link_libraries(BRAINSMushTest BRAINSMushCOMMONLIB)
set_target_properties(BRAINSMushTest PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

##  Test goes here:
set(BRAINSMushTestName BRAINSMushTest_MaskGiven)
add_test(${BRAINSMushTestName} BRAINSMushTest
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/TestData/${BRAINSMushTestName}_MushResult_BrainMask.nii.gz
            ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSMushTest
  --inputFirstVolume ${CMAKE_CURRENT_SOURCE_DIR}/TestData/T1Aligned.nii.gz
  --inputSecondVolume ${CMAKE_CURRENT_SOURCE_DIR}/TestData/T2Aligned.nii.gz
  --inputMaskVolume ${CMAKE_CURRENT_SOURCE_DIR}/TestData/BrainMaskInit.nii.gz
  --outputWeightsFile ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_weights.txt
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_MushTest.nii.gz
  --outputMask ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
)

set(BRAINSMushTestName BRAINSMushTest_NoMaskGiven)
add_test(${BRAINSMushTestName} BRAINSMushTest
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/TestData/${BRAINSMushTestName}_MushResult_BrainMask.nii.gz
            ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSMushTest
  --inputFirstVolume ${CMAKE_CURRENT_SOURCE_DIR}/TestData/T1Aligned.nii.gz
  --inputSecondVolume ${CMAKE_CURRENT_SOURCE_DIR}/TestData/T2Aligned.nii.gz
  --outputWeightsFile ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_weights.txt
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_MushTest.nii.gz
  --outputMask ${CMAKE_CURRENT_BINARY_DIR}/Testing/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
)

### The last section is concerned with installing the binaries and making distributions.
install(TARGETS BRAINSMush
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib/static)

include(InstallRequiredSystemLibraries)


set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "6")
set(CPACK_PACKAGE_VERSION_PATCH "0")

set(CPACK_PACKAGE_NAME "BRAINSMush")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-dev")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CPACK_PACKAGE_NAME}-A program for medical image registration with mutual information metric.")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")

set(CPACK_PACKAGE_DEFAULT_LOCATION "/opt/${CPACK_PACKAGE_NAME}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}")
set(CPACK_PACKAGING_INSTALL_PREFIX "/")
set(CPACK_SET_DESTDIR ON)

set(CPACK_PACKAGE_VENDOR "NAMIC Tool Developed at The University of Iowa")

set(CPACK_SOURCE_GENERATOR "TGZ;TZ")
#set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;\\.*~")
set(CPACK_SOURCE_IGNORE_FILES "")
set(CPACK_PACKAGE_EXECUTABLES "BRAINSMush";"A program for registering medical images with mutual information.")

if(WIN32 AND NOT UNIX)
  set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  set(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
  set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} My Famous Project")
  set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.nitrc.org")
  set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.nitrc.org")
  set(CPACK_NSIS_CONTACT "hans-johnson@uiowa.edu")
  set(CPACK_NSIS_MODIFY_PATH ON)
else(WIN32 AND NOT UNIX)
  set(CPACK_STRIP_FILES OFF)
  set(CPACK_SOURCE_STRIP_FILES OFF)
endif(WIN32 AND NOT UNIX)

include(CPack)
