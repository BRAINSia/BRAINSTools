
#
# list of all tarballs to download
set(DWITarDir ${TestData_DIR}/DWI_TestData_TARBALLS)

set(ALL_DICOM_TARBALLS
${DWITarDir}/DWIConvertGeSignaHdxBMatrixTest.tar.gz
${DWITarDir}/DWIConvertGeSignaHdxBigEndianTest.tar.gz
${DWITarDir}/DWIConvertGeSignaHdxFSLTest.tar.gz
${DWITarDir}/DWIConvertGeSignaHdxTest.tar.gz
${DWITarDir}/DWIConvertGeSignaHdxtBMatrixTest.tar.gz
${DWITarDir}/DWIConvertGeSignaHdxtTest.tar.gz
${DWITarDir}/DWIConvertPhilipsAchieva2BMatrixTest.tar.gz
${DWITarDir}/DWIConvertPhilipsAchieva2Test.tar.gz
${DWITarDir}/DWIConvertPhilipsAchieva3Test.tar.gz
${DWITarDir}/DWIConvertPhilipsAchieva4Test.tar.gz
${DWITarDir}/DWIConvertPhilipsAchieva6Test.tar.gz
${DWITarDir}/DWIConvertPhilipsAchieva7Test.tar.gz
${DWITarDir}/DWIConvertPhilipsAchievaBMatrixTest.tar.gz
${DWITarDir}/DWIConvertPhilipsAchievaBigEndianTest.tar.gz
${DWITarDir}/DWIConvertPhilipsAchievaTest.tar.gz
${DWITarDir}/DWIConvertSiemensTrioSyngo2004A1Test.tar.gz
${DWITarDir}/DWIConvertSiemensTrioSyngo2004A2Test.tar.gz
${DWITarDir}/DWIConvertSiemensTrioTim2Test.tar.gz
${DWITarDir}/DWIConvertSiemensTrioTim3Test.tar.gz
${DWITarDir}/DWIConvertSiemensTrioTimBigEndian1Test.tar.gz
${DWITarDir}/DWIConvertSiemensTrioTimTest.tar.gz
${DWITarDir}/DWIConvertSiemensVerioTest.tar.gz
${DWITarDir}/GeSignaHDx.tar.gz
${DWITarDir}/GeSignaHDxBigEndian.tar.gz
${DWITarDir}/GeSignaHDxt.tar.gz
${DWITarDir}/PhilipsAchieva1.tar.gz
${DWITarDir}/PhilipsAchieva2.tar.gz
${DWITarDir}/PhilipsAchieva3.tar.gz
${DWITarDir}/PhilipsAchieva4.tar.gz
${DWITarDir}/PhilipsAchieva6.tar.gz
${DWITarDir}/PhilipsAchieva7.tar.gz
${DWITarDir}/PhilipsAchievaBigEndian1.tar.gz
${DWITarDir}/SiemensTrio-Syngo2004A-1.tar.gz
${DWITarDir}/SiemensTrio-Syngo2004A-2.tar.gz
${DWITarDir}/SiemensTrioTim1.tar.gz
${DWITarDir}/SiemensTrioTim2.tar.gz
${DWITarDir}/SiemensTrioTim3.tar.gz
${DWITarDir}/SiemensTrioTimBigEndian1.tar.gz
${DWITarDir}/SiemensVerio.tar.gz
)


# accumulator for test names that actually force
# unpacking tars to directories
set(DWITarDownloadTestNames)

foreach( curr_tarball ${ALL_DICOM_TARBALLS} )
  # use base tar name for test name
  get_filename_component(base_filename "${curr_tarball}" NAME)
  set(test_name Uncompress_${base_filename})

  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${test_name}
    WORKING_DIRECTORY ${DWIConvert_BINARY_DIR}
    COMMAND  ${CMAKE_COMMAND} -E tar xzf DATA{${curr_tarball}}
  )
  list(APPEND DWITarDownloadTestNames ${test_name})
endforeach()

#
# make a test depend on the tarball downloads.
# This will make the test TestName depend on ALL the
# test targets responsible for untarring dicom series,
# but it beats having a unique dependency clause for
# each test
function(AddDataDependency TestName)
  foreach(datatest ${DWITarDownloadTestNames})
    set_property(TEST ${TestName} APPEND PROPERTY DEPENDS ${datatest})
  endforeach()
endfunction()


#
# Test CLI executables
#
include_regular_expression("^.*$")

include_directories(${DWIConvert_SOURCE_DIR})
include_directories(${DWIConvert_BINARY_DIR})
include_directories(${ModuleDescriptionParser_SOURCE_DIR})

#
# build compare programs
foreach(CLP DWICompare DWISimpleCompare)
  set ( ${CLP}_SOURCE ${CLP}.cxx)
  generateclp(${CLP}_SOURCE ${CLP}.xml)
  add_executable(${CLP} ${${CLP}_SOURCE})
  target_link_libraries (${CLP} DWIConvertSupportLib BRAINSCommonLib )
endforeach()

## Tool for comparing if BVec and BVal files are the same.
add_executable(BFileCompareTool BFileCompareTool.cpp)
target_link_libraries(BFileCompareTool DWIConvertSupportLib BRAINSCommonLib )

#
# build FSL test program
add_executable(FSLTextFileCompare FSLTextFileCompare.cxx)
target_link_libraries (FSLTextFileCompare DWIConvertSupportLib BRAINSCommonLib )

#
# This program not currently used but hey, might come in
# handy.
add_executable(NrrdToNIfTI NrrdToNIfTI.cxx)
target_link_libraries(NrrdToNIfTI DWIConvertSupportLib BRAINSCommonLib )

#
# This isn't used; we use a CMake Script to handle running the
# tests. The test framework file compare wouldn't work properly
# with vector images
# set (CLP DWIConvert)
# add_executable(${CLP}Test ${CLP}Test.cxx ${DWIConvertTest_SRC})
# target_link_libraries(${CLP}Test DWIConvertSupportLib BRAINSCommonLib )
#
# add_test(NAME DWIConvertHelpTest
#   COMMAND ${SLICER_LAUNCH_COMMAND} $<TARGET_FILE:DWIConvertTest> DWIConvertTest --help
#   )

# dump results in current binary dir.
set(TstOutput ${CMAKE_CURRENT_BINARY_DIR})

set(DWIBASELINE_DIR ${TestData_DIR}/DWI_TestData_OUTPUTS)

#
# pretty useless test except for slight increase in coverage
# and verify the programs got built
add_test(NAME DWIConvertDWICompareHelpTest
  COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:DWICompare>
  --help
  )

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertGeSignaHdxTest COMMAND ${CMAKE_COMMAND}
  -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/GeSignaHDx.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/GeSignaHDx
  -D TEST_TEMP_OUTPUT=${TstOutput}/GeSignaHDxTest.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertGeSignaHdxTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertGeSignaHdxFSLTest COMMAND ${CMAKE_COMMAND}
  -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWISimpleCompare>
  -D TEST_TEXT_COMPARE=$<TARGET_FILE:FSLTextFileCompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/GeSignaHDx.nii.gz}
  -D TEST_BASELINE_BVEC=DATA{${DWIBASELINE_DIR}/GeSignaHDx.bvec}
  -D TEST_BASELINE_BVAL=DATA{${DWIBASELINE_DIR}/GeSignaHDx.bval}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/GeSignaHDx
  -D TEST_TEMP_OUTPUT=${TstOutput}/GeSignaHDxTest.nii.gz
  -D TEST_FSL_FLAG=ON
  -D TEST_TEMP_BVEC=${TstOutput}/GeSignaHDxTest.bvec
  -D TEST_TEMP_BVAL=${TstOutput}/GeSignaHDxTest.bval
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertGeSignaHdxFSLTest)

add_test(NAME FSLToNrrdTest
  COMMAND ${CMAKE_COMMAND}
  -DVEC_FILE=${TstOutput}/GeSignaHDxTest.bvec
  -DVAL_FILE=${TstOutput}/GeSignaHDxTest.bval
  -DNII_FILE=${TstOutput}/GeSignaHDxTest.nii.gz
  -DNRRD_FILE=${TstOutput}/FSLToNrrdTest.nrrd
  -DNRRD_COMPARE_FILE=${TstOutput}/GeSignaHDxTest.nrrd
  -DFSL_TO_NRRD=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -P ${CMAKE_CURRENT_LIST_DIR}/FSLToNrrdTest.cmake
  )

set_property(TEST FSLToNrrdTest
  APPEND PROPERTY DEPENDS DWIConvertGeSignaHdxFSLTest
                          DWIConvertGeSignaHdxTest)

add_test(NAME NrrdToFSLTest
  COMMAND ${CMAKE_COMMAND}
  -DNRRD_FILE=${TstOutput}/FSLToNrrdTest.nrrd
  -DVEC_COMPARE_FILE=${TstOutput}/GeSignaHDxTest.bvec
  -DVAL_COMPARE_FILE=${TstOutput}/GeSignaHDxTest.bval
  -DNII_COMPARE_FILE=${TstOutput}/GeSignaHDxTest.nii.gz
  -DVEC_FILE=${TstOutput}/NrrdToFSLTest.bvec
  -DVAL_FILE=${TstOutput}/NrrdToFSLTest.bval
  -DNII_FILE=${TstOutput}/NrrdToFSLTest.nii.gz
  -DNRRD_TO_FSL=$<TARGET_FILE:DWIConvert>
  -DTEST_COMPARE_PROGRAM=$<TARGET_FILE:DWISimpleCompare>
  -DTEXT_COMPARE_PROGRAM=$<TARGET_FILE:FSLTextFileCompare>
  -P ${CMAKE_CURRENT_LIST_DIR}/NrrdToFSLTest.cmake
  )

set_property(TEST NrrdToFSLTest
  APPEND PROPERTY DEPENDS FSLToNrrdTest)

# this tests the case where DWIConvert should generate the bval/bvec
# names based on the Output NIfTI volume
add_test(NAME NrrdToFSLTest2
  COMMAND ${CMAKE_COMMAND} --trace
  -DNRRD_FILE=${TstOutput}/FSLToNrrdTest.nrrd
  -DVEC_COMPARE_FILE=${TstOutput}/GeSignaHDxTest.bvec
  -DVAL_COMPARE_FILE=${TstOutput}/GeSignaHDxTest.bval
  -DNII_COMPARE_FILE=${TstOutput}/GeSignaHDxTest.nii.gz
  -DNII_FILE=${TstOutput}/NrrdToFSLTest2.nii.gz
  -DNRRD_TO_FSL=$<TARGET_FILE:DWIConvert>
  -DTEST_COMPARE_PROGRAM=$<TARGET_FILE:DWISimpleCompare>
  -DTEXT_COMPARE_PROGRAM=$<TARGET_FILE:FSLTextFileCompare>
  -P ${CMAKE_CURRENT_LIST_DIR}/NrrdToFSLTest.cmake
  )

set_property(TEST NrrdToFSLTest2
  APPEND PROPERTY DEPENDS FSLToNrrdTest)


ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertGeSignaHdxBigEndianTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/GeSignaHDx.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/GeSignaHDxBigEndian
  -D TEST_TEMP_OUTPUT=${TstOutput}/GeSignaHDxBigEndianTest.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertGeSignaHdxBigEndianTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertGeSignaHdxBMatrixTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/GeSignaHDx.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/GeSignaHDx
  -D TEST_TEMP_OUTPUT=${TstOutput}/GeSignaHDxTest.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertGeSignaHdxBMatrixTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertGeSignaHdxtTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/GeSignaHDxt.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/GeSignaHDxt
  -D TEST_TEMP_OUTPUT=${TstOutput}/GeSignaHDxtTest.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertGeSignaHdxtTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertGeSignaHdxtBMatrixTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/GeSignaHDxt.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/GeSignaHDxt
  -D TEST_TEMP_OUTPUT=${TstOutput}/GeSignaHDxtTest.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertGeSignaHdxtBMatrixTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchievaTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva1.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva1
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva1Test.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchievaTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchievaBigEndianTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva1.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchievaBigEndian1
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchievaBigEndian1Test.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchievaBigEndianTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchievaBMatrixTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva1.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva1
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva1Test.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchievaBMatrixTest)
set_tests_properties(DWIConvertPhilipsAchievaBMatrixTest PROPERTIES DEPENDS DWIConvertPhilipsAchievaTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchieva2Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva2.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva2
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva2Test.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchieva2Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchieva2BMatrixTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva2.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva2
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva2Test.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchieva2BMatrixTest)
set_tests_properties(DWIConvertPhilipsAchieva2BMatrixTest PROPERTIES DEPENDS DWIConvertPhilipsAchieva2Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchieva3Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva3.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva3
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva3Test.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchieva3Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchieva4Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva4.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva4
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva4Test.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchieva4Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchieva6Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva6.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva6
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva6Test.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchieva6Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertPhilipsAchieva7Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/PhilipsAchieva7.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/PhilipsAchieva7
  -D TEST_TEMP_OUTPUT=${TstOutput}/PhilipsAchieva7Test.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertPhilipsAchieva7Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertSiemensTrioTimTest
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/SiemensTrioTim1.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/SiemensTrioTim1
  -D TEST_TEMP_OUTPUT=${TstOutput}/SiemensTrioTim1Test.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertSiemensTrioTimTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertSiemensTrioTimBigEndian1Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/SiemensTrioTim1.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/SiemensTrioTimBigEndian1
  -D TEST_TEMP_OUTPUT=${TstOutput}/SiemensTrioTimBigEndian1Test.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertSiemensTrioTimBigEndian1Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertSiemensTrioTim2Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/SiemensTrioTim2.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/SiemensTrioTim2
  -D TEST_TEMP_OUTPUT=${TstOutput}/SiemensTrioTim2Test.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertSiemensTrioTim2Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertSiemensTrioTim3Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/SiemensTrioTim3.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/SiemensTrioTim3
  -D TEST_TEMP_OUTPUT=${TstOutput}/SiemensTrioTim3Test.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertSiemensTrioTim3Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertSiemensVerioTest
  COMMAND ${CMAKE_COMMAND}
  -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/SiemensVerio.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/SiemensVerio
  -D TEST_TEMP_OUTPUT=${TstOutput}/SiemensVerio.nrrd
  -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertSiemensVerioTest)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertSiemensTrio_Syngo2004A_1Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-1.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/SiemensTrio-Syngo2004A-1
  -D TEST_TEMP_OUTPUT=${TstOutput}/SiemensTrio-Syngo2004A-1.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertSiemensTrio_Syngo2004A_1Test)

ExternalData_Add_Test( ${PROJECT_NAME}FetchData NAME DWIConvertSiemensTrio_Syngo2004A_2Test
  COMMAND ${CMAKE_COMMAND} -D TEST_PROGRAM=$<TARGET_FILE:DWIConvert>
  -D TEST_COMPARE_PROGRAM=$<TARGET_FILE:DWICompare>
  -D TEST_BASELINE=DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-2.nrrd}
  -D TEST_INPUT=${DWIConvert_BINARY_DIR}/SiemensTrio-Syngo2004A-2
  -D TEST_TEMP_OUTPUT=${TstOutput}/SiemensTrio-Syngo2004A-2.nrrd
  -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )
AddDataDependency(DWIConvertSiemensTrio_Syngo2004A_2Test)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
