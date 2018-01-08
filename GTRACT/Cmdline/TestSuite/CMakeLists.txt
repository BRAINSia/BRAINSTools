
#add_executable( gtractCoRegAnatomyTests
#    ${GTRACT_SOURCE_DIR}/Cmdline/gtractCoRegAnatomyPrimary.cxx gtractCoRegAnatomyTests.cxx )
#target_link_libraries( gtractCoRegAnatomyTests GTRACTCommon BRAINSCommonLib DWIConvertSupportLib )
#set_target_properties(gtractCoRegAnatomyTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${GTRACT_BINARY_DIR})
#set(GTRACT_CoRegAnatomy_TESTS ${GTRACT_BINARY_DIR}/gtractCoRegAnatomyTests )
#
#
#set(GTRACT_TEST_OUTPUT_DIR "${GTRACT_BINARY_DIR}/Testing")
#make_directory(${GTRACT_TEST_OUTPUT_DIR})
#
#set(GTRACT_TEST_INPUT_DIR "${GTRACT_SOURCE_DIR}/Data")
#
#message(warning "${CMAKE_BUNDLE_OUTPUT_DIRECTORY}")
#
#add_test(NAME gtractCoRegAnatomy_Rigid_Test COMMAND ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoRegAnatomy_TESTS>
#   gtractCoRegAnatomyTest
#     --inputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
#     --vectorIndex 0
#     --inputAnatomicalVolume ${GTRACT_TEST_INPUT_DIR}/t1_002.nhdr
#     --outputTransformName ${GTRACT_TEST_OUTPUT_DIR}/rigid_s.${XFRM_EXT}
#     --numberOfSamples 500000
#     --transformType Rigid
#)
#
#add_test(NAME gtractCoRegAnatomy_Bspline_Test COMMAND ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoRegAnatomy_TESTS>
#   gtractCoRegAnatomyTest
#     --inputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
#     --vectorIndex 0
#     --inputAnatomicalVolume ${GTRACT_TEST_INPUT_DIR}/t1_002.nhdr
#     --outputTransformName ${GTRACT_TEST_OUTPUT_DIR}/bspline_s.${XFRM_EXT}
#     --numberOfIterations 500
#     --transformType Bspline
#     --gradientTolerance 5e-6
#     --inputRigidTransform ${GTRACT_TEST_OUTPUT_DIR}/rigid_s.${XFRM_EXT}
#     --gridSize 18,6,18
#     --spatialScale 10
#)

set(DWIBASELINE_DIR ${TestData_DIR}/DWI_TestData_OUTPUTS)

## Test for gtractCoregBvalues
add_executable( gtractCoregBvaluesTests gtractCoregBvaluesTests.cxx )
target_link_libraries( gtractCoregBvaluesTests GTRACTCommon BRAINSCommonLib DWIConvertSupportLib)
set_target_properties(gtractCoregBvaluesTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

## TODO: gtractCoregBvalues returns segfault after the registration step.
##       This test is commented out until the bug is fixed.
##       The use of this program seems to be deprecated, since DTIPrep is
##       used for coregistration and eddy current correction.
#set(GTRACTTestName gtractCoregBvalues_Test)
#ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
#  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractCoregBvaluesTests>
#    --compare DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-1.nrrd}
#              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
#    --compareNumberOfPixelsTolerance 1000
#    --compareIntensityTolerance 20
# gtractCoregBvaluesTest
#  --fixedVolume DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-1.nrrd}
#  --movingVolume DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-1.nrrd}
#  --fixedVolumeIndex 0
#  --debugLevel 0
#  --eddyCurrentCorrection
#  --relaxationFactor 0.25
#  --maximumStepSize 0.1
#  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
#  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.${XFRM_EXT}
#)

## Test for gtractAverageBvalues
add_executable( gtractAverageBvaluesTests gtractAverageBvaluesTests.cxx )
target_link_libraries( gtractAverageBvaluesTests BRAINSCommonLib GTRACTCommon DWIConvertSupportLib)
set_target_properties(gtractAverageBvaluesTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(GTRACTTestName GTRACTTest_gtractAverageBvalues)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractAverageBvaluesTests>
  --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
            ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractAverageBvaluesTest
    --inputVolume DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-1.nrrd}
    --averageB0only
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
)

## Test for gtractConcatDwi
add_executable( gtractConcatDwiTests gtractConcatDwiTests.cxx )
target_link_libraries( gtractConcatDwiTests BRAINSCommonLib GTRACTCommon DWIConvertSupportLib)
set_target_properties(gtractConcatDwiTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(GTRACTTestName GTRACTTest_gtractConcatDwi)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractConcatDwiTests>
  --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
            ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractConcatDwiTest
    --inputVolume DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-1.nrrd},DATA{${DWIBASELINE_DIR}/SiemensTrio-Syngo2004A-2.nrrd}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
)

# add the directory containing DWI Dicom NiFTI files
set(DWITestFileDir ${PROJECT_BINARY_DIR}/../DWIConvert)
#set(DWITestFileDir /Users/hxie1/Projects/BrainsToolsBuild/DWIConvert)

# in order to make sure the origin are same of concating images, we use same images to concatenate

# new added test case for GTRACT to support NifTI
# For NifTI, need to find suitible NifTI file for test
#add_test(NAME ${GTRACTTestName}_Concat_NifTI
#        COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractConcatDwiTests>
#        gtractConcatDwiTest
#        --inputVolume ${DWITestFileDir}/DWIConvertMultiShellFSITest_Recover.nii,${DWITestFileDir}/DWIConvertMultiShellFSITest_Recover.nii
#        --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.nii.test.nrrd
#        )

# new added test case for GTRACT to support Dicom
add_test(NAME ${GTRACTTestName}_Concat_Dicom
        COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractConcatDwiTests>
        gtractConcatDwiTest
        --inputVolume ${DWITestFileDir}/PhilipsAchieva1,${DWITestFileDir}/PhilipsAchieva1
        --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.dicom.test.nrrd
        )

# new added test case for GTRACT to support multi Dicom
add_test(NAME ${GTRACTTestName}_Concat_multi_Dicom
        COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractConcatDwiTests>
        gtractConcatDwiTest
        --inputVolume ${DWITestFileDir}/PhilipsAchieva1,${DWITestFileDir}/PhilipsAchieva1,${DWITestFileDir}/PhilipsAchieva1
        --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.multidicom.test.nrrd
        )


## A set of tests for gtractResampleDWIInPlace
add_executable( gtractResampleDWIInPlaceTests gtractResampleDWIInPlaceTests.cxx )
target_link_libraries( gtractResampleDWIInPlaceTests  BRAINSCommonLib GTRACTCommon DWIConvertSupportLib)
set_target_properties(gtractResampleDWIInPlaceTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(GTRACTTestName GTRACTTest_gtractResampleDWIInPlace)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractResampleDWIInPlaceTests>
    --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractResampleDWIInPlaceTest
    --inputVolume DATA{${TestData_DIR}/DWI_test1.nrrd}
    --inputTransform DATA{${TestData_DIR}/${GTRACTTestName}_rigid_transform.mat}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compareTwoCSVFiles.py.in ${CMAKE_CURRENT_BINARY_DIR}/compareTwoCSVFiles.py @ONLY IMMEDIATE)

## The following set of tests verify that gtractResampleDWIInPlace
# changes the gradient directions table in a proper way when the
# input DWI image goes under different rotations.
# 5 tests are designed. The first three tests apply 3 different rotations around z axis (-15, 15, 45);
# then, two more tests are designed by 30 degree rotation around each of x, y directions.
# For each regression test the rotated DWI header file (measurement frame and diffusion
# gradient directions) is written to a simple CSV file that is compared with corresponding
# baseline csv file.
# The reference metaData are approved by looking into corresponding color tensor images
# created by ditestim.

set(GTRACTTestName GTRACTTest_gtractResampleDWIInPlace_Rot15z)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractResampleDWIInPlaceTests>
    --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractResampleDWIInPlaceTest
    --inputVolume DATA{${TestData_DIR}/DWI_test1.nrrd}
    --inputTransform DATA{${TestData_DIR}/GTRACTTest_rigidTx_Rot15z.h5}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --writeOutputMetaData ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}_metaData.test.csv
)

set(PARENT_TEST GTRACTTest_gtractResampleDWIInPlace_Rot15z)
set(GTRACTTestName chk_${PARENT_TEST})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} python compareTwoCSVFiles.py
    --csvReferenceFile DATA{${TestData_DIR}/${GTRACTTestName}_metaData.standard.csv}
    --csvInputFile ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_TEST}_metaData.test.csv
)
set_property(TEST ${GTRACTTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(GTRACTTestName GTRACTTest_gtractResampleDWIInPlace_Rot-15z)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractResampleDWIInPlaceTests>
    --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractResampleDWIInPlaceTest
    --inputVolume DATA{${TestData_DIR}/DWI_test1.nrrd}
    --inputTransform DATA{${TestData_DIR}/GTRACTTest_rigidTx_Rot-15z.h5}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --writeOutputMetaData ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}_metaData.test.csv
)

set(PARENT_TEST GTRACTTest_gtractResampleDWIInPlace_Rot-15z)
set(GTRACTTestName chk_${PARENT_TEST})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} python compareTwoCSVFiles.py
    --csvReferenceFile DATA{${TestData_DIR}/${GTRACTTestName}_metaData.standard.csv}
    --csvInputFile ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_TEST}_metaData.test.csv
)
set_property(TEST ${GTRACTTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(GTRACTTestName GTRACTTest_gtractResampleDWIInPlace_Rot45z)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractResampleDWIInPlaceTests>
    --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractResampleDWIInPlaceTest
    --inputVolume DATA{${TestData_DIR}/DWI_test1.nrrd}
    --inputTransform DATA{${TestData_DIR}/GTRACTTest_rigidTx_Rot45z.h5}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --writeOutputMetaData ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}_metaData.test.csv
)

set(PARENT_TEST GTRACTTest_gtractResampleDWIInPlace_Rot45z)
set(GTRACTTestName chk_${PARENT_TEST})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} python compareTwoCSVFiles.py
    --csvReferenceFile DATA{${TestData_DIR}/${GTRACTTestName}_metaData.standard.csv}
    --csvInputFile ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_TEST}_metaData.test.csv
)
set_property(TEST ${GTRACTTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(GTRACTTestName GTRACTTest_gtractResampleDWIInPlace_Rot30x)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractResampleDWIInPlaceTests>
    --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractResampleDWIInPlaceTest
    --inputVolume DATA{${TestData_DIR}/DWI_test1.nrrd}
    --inputTransform DATA{${TestData_DIR}/GTRACTTest_rigidTx_Rot30x.h5}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --writeOutputMetaData ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}_metaData.test.csv
)

set(PARENT_TEST GTRACTTest_gtractResampleDWIInPlace_Rot30x)
set(GTRACTTestName chk_${PARENT_TEST})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} python compareTwoCSVFiles.py
    --csvReferenceFile DATA{${TestData_DIR}/${GTRACTTestName}_metaData.standard.csv}
    --csvInputFile ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_TEST}_metaData.test.csv
)
set_property(TEST ${GTRACTTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(GTRACTTestName GTRACTTest_gtractResampleDWIInPlace_Rot30y)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractResampleDWIInPlaceTests>
    --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractResampleDWIInPlaceTest
    --inputVolume DATA{${TestData_DIR}/DWI_test1.nrrd}
    --inputTransform DATA{${TestData_DIR}/GTRACTTest_rigidTx_Rot30y.h5}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --writeOutputMetaData ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}_metaData.test.csv
)

set(PARENT_TEST GTRACTTest_gtractResampleDWIInPlace_Rot30y)
set(GTRACTTestName chk_${PARENT_TEST})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} python compareTwoCSVFiles.py
    --csvReferenceFile DATA{${TestData_DIR}/${GTRACTTestName}_metaData.standard.csv}
    --csvInputFile ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_TEST}_metaData.test.csv
)
set_property(TEST ${GTRACTTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})


ExternalData_Add_Target( ${PROJECT_NAME}FetchData )