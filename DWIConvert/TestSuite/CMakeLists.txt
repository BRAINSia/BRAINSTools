#
# Test CLI executables
#
include_regular_expression("^.*$")

include_directories(${DWIConvert_SOURCE_DIR})
include_directories(${DWIConvert_BINARY_DIR})
include_directories(${ModuleDescriptionParser_SOURCE_DIR})


foreach(CLP DWICompare DWISimpleCompare)
  set ( ${CLP}_SOURCE ${CLP}.cxx)
  generateclp(${CLP}_SOURCE ${CLP}.xml)
  add_executable(${CLP} ${${CLP}_SOURCE})
  target_link_libraries (${CLP} DWIConvertSupportLib )
  set(${CLP}EXE ${DWIConvert_BINARY_DIR}/TestSuite/${CLP})
endforeach()

add_executable(FSLTextFileCompare FSLTextFileCompare.cxx)
target_link_libraries (FSLTextFileCompare DWIConvertSupportLib )
set(FSLTextFileCompareEXE
  ${DWIConvert_BINARY_DIR}/TestSuite/FSLTextFileCompare)

add_executable(NrrdToNIfTI NrrdToNIfTI.cxx)
target_link_libraries(NrrdToNIfTI DWIConvertSupportLib)

set (CLP DWIConvert)
add_executable(${CLP}Test ${CLP}Test.cxx ${DWIConvertTest_SRC})
add_dependencies(${CLP}Test ${CLP})
target_link_libraries(${CLP}Test DWIConvertSupportLib )

set (DWIConvertEXE ${Slicer_LAUNCH_COMMAND} ${DWIConvert_BINARY_DIR}/${CLP} )
set (DWIConvert_TESTS ${Slicer_LAUNCH_COMMAND} ${DWIConvert_BINARY_DIR}/TestSuite/${CLP}Test )

if(NOT SETIFEMPTY)
macro(SETIFEMPTY)
  set(KEY ${ARGV0})
  set(VALUE ${ARGV1})
  if(NOT ${KEY})
    set(${ARGV})
  endif(NOT ${KEY})
endmacro(SETIFEMPTY KEY VALUE)
endif(NOT SETIFEMPTY)

set(TEMP ${CMAKE_CURRENT_BINARY_DIR})

include(${CMAKE_CURRENT_LIST_DIR}/MIDAS.cmake)
SETIFEMPTY(MIDAS_REST_URL "http://midas.kitware.com/api/rest" CACHE STRING "The MIDAS server where testing data resides")

set(MIDAS_KEY_DIR "${DWIConvert_SOURCE_DIR}/TestSuite/keys")

add_test(DWIConvertHelpTest ${DWIConvert_TESTS}
    DWIConvertTest
    --help
  )

add_test(DWIConvertDWICompareHelpTest ${Slicer_LAUNCH_COMMAND} ${DWICompareEXE}
    --help
  )

midas_add_test(NAME DWIConvertGeSignaHdxTest COMMAND ${CMAKE_COMMAND}
        -D TEST_PROGRAM=${DWIConvertEXE}
        -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
        -D TEST_BASELINE=MIDAS{GeSignaHDx.nrrd.md5}
        -D TEST_INPUT=MIDAS_TGZ{GeSignaHDx.tar.gz.md5}
        -D TEST_TEMP_OUTPUT=${TEMP}/GeSignaHDxTest.nrrd
        -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
        )

midas_add_test(NAME DWIConvertGeSignaHdxFSLTest COMMAND ${CMAKE_COMMAND}
        -D TEST_PROGRAM=${DWIConvertEXE}
        -D TEST_COMPARE_PROGRAM=${DWISimpleCompareEXE}
        -D TEST_BASELINE=MIDAS{GeSignaHDx.nii.gz.md5}
        -D TEST_BASELINE_BVEC=MIDAS{GeSignaHDx.bvec.md5}
        -D TEST_BASELINE_BVAL=MIDAS{GeSignaHDx.bval.md5}
        -D TEST_INPUT=MIDAS_TGZ{GeSignaHDx.tar.gz.md5}
        -D TEST_TEMP_OUTPUT=${TEMP}/GeSignaHDxTest.nii.gz
        -D TEST_FSL_FLAG=ON
        -D TEST_TEMP_BVEC=${TEMP}/GeSignaHDxTest.bvec
        -D TEST_TEMP_BVAL=${TEMP}/GeSignaHDxTest.bval
        -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
        )

add_test(NAME FSLToNrrdTest
  COMMAND ${CMAKE_COMMAND}
  -DVEC_FILE=${TEMP}/GeSignaHDxTest.bvec
  -DVAL_FILE=${TEMP}/GeSignaHDxTest.bval
  -DNII_FILE=${TEMP}/GeSignaHDxTest.nii.gz
  -DNRRD_FILE=${TEMP}/FSLToNrrdTest.nrrd
  -DNRRD_COMPARE_FILE=${TEMP}/GeSignaHDxTest.nrrd
  -DFSL_TO_NRRD=${DWIConvert_BINARY_DIR}/DWIConvert
  -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
  -P ${CMAKE_CURRENT_LIST_DIR}/FSLToNrrdTest.cmake
)

set_property(TEST FSLToNrrdTest
  APPEND PROPERTY DEPENDS DWIConvertGeSignaHdxFSLTest)

add_test(NAME NrrdToFSLTest
  COMMAND ${CMAKE_COMMAND}
  -DNRRD_FILE=${TEMP}/FSLToNrrdTest.nrrd
  -DVEC_COMPARE_FILE=${TEMP}/GeSignaHDxTest.bvec
  -DVAL_COMPARE_FILE=${TEMP}/GeSignaHDxTest.bval
  -DNII_COMPARE_FILE=${TEMP}/GeSignaHDxTest.nii.gz
  -DVEC_FILE=${TEMP}/NrrdToFSLTest.bvec
  -DVAL_FILE=${TEMP}/NrrdToFSLTest.bval
  -DNII_FILE=${TEMP}/NrrdToFSLTest.nii.gz
  -DNRRD_TO_FSL=${DWIConvert_BINARY_DIR}/DWIConvert
  -DTEST_COMPARE_PROGRAM=${DWISimpleCompareEXE}
  -DTEXT_COMPARE_PROGRAM=${FSLTextFileCompareEXE}
  -P ${CMAKE_CURRENT_LIST_DIR}/NrrdToFSLTest.cmake
)

set_property(TEST NrrdToFSLTest
  APPEND PROPERTY DEPENDS FSLToNrrdTest)

midas_add_test(NAME DWIConvertGeSignaHdxBigEndianTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{GeSignaHDx.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{GeSignaHDxBigEndian.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/GeSignaHDxBigEndianTest.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertGeSignaHdxBMatrixTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{GeSignaHDx.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{GeSignaHDx.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/GeSignaHDxTest.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertGeSignaHdxtTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{GeSignaHDxt.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{GeSignaHDxt.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/GeSignaHDxtTest.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertGeSignaHdxtBMatrixTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{GeSignaHDxt.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{GeSignaHDxt.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/GeSignaHDxtTest.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchievaTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva1.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva1.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva1Test.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchievaBigEndianTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva1.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchievaBigEndian1.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchievaBigEndian1Test.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchievaBMatrixTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva1.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva1.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva1Test.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(DWIConvertPhilipsAchieva2Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva2.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva2.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva2Test.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchieva2BMatrixTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva2.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva2.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva2Test.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchieva3Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva3.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva3.tar.gz.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva3Test.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchieva4Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva4.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva4.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva4Test.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchieva6Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva6.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva6.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva6Test.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertPhilipsAchieva7Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{PhilipsAchieva7.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{PhilipsAchieva7.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/PhilipsAchieva7Test.nrrd
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertSiemensTrioTimTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{SiemensTrioTim1.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{SiemensTrioTim1.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/SiemensTrioTim1Test.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertSiemensTrioTimBigEndian1Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{SiemensTrioTim1.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{SiemensTrioTimBigEndian1.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/SiemensTrioTimBigEndian1Test.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertSiemensTrioTim2Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{SiemensTrioTim2.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{SiemensTrioTim2.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/SiemensTrioTim2Test.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertSiemensTrioTim3Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{SiemensTrioTim3.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{SiemensTrioTim3.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/SiemensTrioTim3Test.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertSiemensVerioTest COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{SiemensVerio.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{SiemensVerio.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/SiemensVerio.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertSiemensTrio_Syngo2004A_1Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{SiemensTrio-Syngo2004A-1.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{SiemensTrio-Syngo2004A-1.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/SiemensTrio-Syngo2004A-1.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

midas_add_test(NAME DWIConvertSiemensTrio_Syngo2004A_2Test COMMAND ${CMAKE_COMMAND}
  ${CMAKE_COMMAND} -D TEST_PROGRAM=${DWIConvertEXE}
                   -D TEST_COMPARE_PROGRAM=${DWICompareEXE}
                   -D TEST_BASELINE=MIDAS{SiemensTrio-Syngo2004A-2.nrrd.md5}
                   -D TEST_INPUT=MIDAS_TGZ{SiemensTrio-Syngo2004A-2.md5}
                   -D TEST_TEMP_OUTPUT=${TEMP}/SiemensTrio-Syngo2004A-2.nrrd
                   -D TEST_PROGRAM_ARGS=--useBMatrixGradientDirections
                   -P ${CMAKE_CURRENT_LIST_DIR}/DicomToNrrdDWICompareTest.cmake
  )

