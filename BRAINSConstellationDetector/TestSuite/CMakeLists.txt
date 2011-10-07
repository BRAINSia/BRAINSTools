## Set testing environment
##
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/obliq_setup.txt.in
  ${CMAKE_CURRENT_BINARY_DIR}/obliq_setup.txt @ONLY IMMEDIATE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/T1.mdl ${CMAKE_CURRENT_BINARY_DIR}/T1.mdl)

set(ALL_TEST_PROGS
  BRAINSAlignMSP
  BRAINSConstellationDetector
  BRAINSEyeDetector
  BRAINSClipInferior
  BRAINSConstellationModeler
  BRAINSTrimForegroundInDirection
  )


foreach(testprog ${ALL_TEST_PROGS})
  MakeTestDriverFromSEMTool(${testprog} ${testprog}Test.cxx)
endforeach()


## Test BRAINSConstellationModeler
##

if ( 0 ) ## HACK:  Just silencing failing tests until they can be proper addressed
  set(BRAINSConstellationModelerTestName BRAINSConstellationModelerTest_T1)
  if(ConstellationModeler_TEST_causes_major_delay)
    ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationModelerTestName}
      COMMAND BRAINSConstellationModelerTest
      BRAINSConstellationModelerTest
      --inputTrainingList ${CMAKE_CURRENT_BINARY_DIR}/obliq_setup.txt
      --outputModel ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationModelerTestName}.mdl
      )
  endif(ConstellationModeler_TEST_causes_major_delay)
endif(0)

## Test BRAINSEyeDetector
##
set(BRAINSEyeDetectorTestName BRAINSEyeDetectorTest_T1)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSEyeDetectorTestName}
  COMMAND $<TARGET_FILE:BRAINSEyeDetectorTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSEyeDetectorTestName}_standard.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSEyeDetectorTestName}_aligned.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSEyeDetectorTest
  --inputVolume DATA{${TestData_DIR}/${BRAINSEyeDetectorTestName}_standard.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSEyeDetectorTestName}_aligned.nii.gz
  --debugDir ${CMAKE_CURRENT_BINARY_DIR})

## Test BRAINSAlignMSP:
##
set(BRAINSAlignMSPTestName BRAINSAlignMSPTest_T1)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSAlignMSPTestName}
  COMMAND $<TARGET_FILE:BRAINSAlignMSPTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSAlignMSPTestName}_standard.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSAlignMSPTestName}_aligned.nii.gz
  --compareIntensityTolerance 50
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 500
  BRAINSAlignMSPTest
  --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
  --OutputresampleMSP ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSAlignMSPTestName}_aligned.nii.gz
  --mspQualityLevel 3
  )

## Test BRAINSConstellationDetector
##


if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_acLowerBound)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --acLowerBound 100)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_backgroundValue)

  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --BackgroundFillValue BIGNEG)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_cutOutHeadInOutputResampledVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --cutOutHeadInOutputVolume)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_cutOutHeadInOutputVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --cutOutHeadInOutputVolume)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_interpolationMode)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --interpolationMode BSpline)


  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_mspLevel)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --mspQualityLevel 3)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_otsuThreshold)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --otsuPercentileThreshold 0.05)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputUntransformedClippedVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputUntransformedClippedVolume
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --acLowerBound 100
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5})

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rac)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --rac 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rescaleIntensity)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --rescaleIntensities)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rescaleIntensityRange)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --rescaleIntensitiesOutputRange 20,2000)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rmpj)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --rmpj 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rpc)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --rpc 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rVN4)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --rVN4 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rVN4-rpc-rac-rmpj)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --rmpj 200
    --rVN4 200
    --rac 200
    --rpc 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_T1_outputResampledVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5})

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_T1_outputVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5})

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_trimRescaledIntensities)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
    --trimRescaledIntensities 2)

endif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_forceHoughReportFailure)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/EMSP.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
  --inputVolume DATA{${TestData_DIR}/2.nii.gz}
  --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
  --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
  --resultsDir ${CMAKE_CURRENT_BINARY_DIR}/
  --forceHoughEyeDetectorReportFailure)
set_property(TEST ${BRAINSConstellationDetectorTestName}
            PROPERTY WILL_FAIL TRUE)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_inputLandmarksEMSP)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
  --inputVolume DATA{${TestData_DIR}/2.nii.gz}
  --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
  --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
  --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
  --inputLandmarksEMSP DATA{${TestData_DIR}/EMSP.fcsv})



set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_inputVolumeError)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>

        BRAINSConstellationDetectorTest)
set_property(TEST ${BRAINSConstellationDetectorTestName}
            PROPERTY WILL_FAIL TRUE)

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_LLSModelError)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>

        BRAINSConstellationDetectorTest
        --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
        --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz)
set_property(TEST ${BRAINSConstellationDetectorTestName}
            PROPERTY WILL_FAIL TRUE)

## Test BRAINSTrimForegroundInDirection:
##
set(BRAINSTrimForegroundInDirectionTestName BRAINSTrimForegroundInDirectionTest_T1)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSTrimForegroundInDirectionTestName}
  COMMAND $<TARGET_FILE:BRAINSTrimForegroundInDirectionTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSTrimForegroundInDirectionTestName}_standard.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSTrimForegroundInDirectionTestName}_aligned.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0

  BRAINSTrimForegroundInDirectionTest
  --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSTrimForegroundInDirectionTestName}_aligned.nii.gz
  --headSizeLimit 2600.0)


## Test BRAINSClipInferior
##
# add_executable(BRAINSClipInferiorTest BRAINSClipInferiorTest.cxx)
# target_link_libraries(BRAINSClipInferiorTest BRAINSClipInferiorCOMMONLIB)


## Test Versor
##
set (VersorTester_source VersorTester.cxx)
add_executable(VersorTester ${VersorTester_source})
target_link_libraries(VersorTester ${ITK_LIBRARIES})

install(TARGETS VersorTester
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static)


## Test landmarksConstellationTrainingDefinitionIO
##
add_executable(TestlandmarksConstellationTrainingDefinitionIO TestlandmarksConstellationTrainingDefinitionIO.cxx)


#
# Test Landmark files
#
include_regular_expression("^.*$")

enable_testing()

if(NOT ITK_FOUND)
  find_package(ITK REQUIRED)
  include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

include_directories(${BRAINSConstellationDetectorTool_SOURCE_DIR})
include_directories(${BRAINSConstellationDetectorTool_BINARY_DIR})

find_package(GenerateCLP NO_MODULE REQUIRED)
include(${GenerateCLP_USE_FILE})

set (CLP LandmarksCompare)
set ( ${CLP}_SOURCE ${CLP}.cxx)
generateclp(${CLP}_SOURCE ${CLP}.xml)
add_library(${CLP}Lib SHARED ${${CLP}_SOURCE})
add_executable(${CLP} ${${CLP}_SOURCE})
target_link_libraries(${CLP} ${ITK_LIBRARIES})
set_target_properties(${CLP} PROPERTIES
                        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(LandmarksCompareEXE ${CMAKE_CURRENT_BINARY_DIR}/${CLP})
set (CLP BRAINSConstellationDetector)
set (BRAINSConstellationDetectorEXE ${LAUNCH_EXE} ${BRAINSConstellationDetectorTool_BINARY_DIR}/../../bin/${CLP} )
set (BRAINSConstellationDetector_TESTS ${LAUNCH_EXE} ${CMAKE_CURRENT_BINARY_DIR}/${CLP}TestDriver )


if(0)
#This test is run without any output to check the warning codes
ExternalData_add_test(${PROJECT_NAME}FetchData NAME BRAINSConstellationDetectorTestWarning1
  COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
  BRAINSConstellationDetectorTest
  --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
  --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
  --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
  --inputLandmarksPaired DATA{${TestData_DIR}/BRAINSConstellationDetectorTest_outputLandmarksInInputSpace_standard.fcsv
  )
set_tests_properties(BRAINSConstellationDetectorTestWarning1 PROPERTIES WILL_FAIL ON)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME BRAINSConstellationDetectorTestWarning2
  COMMAND $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
  BRAINSConstellationDetectorTest
  --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
  --LLSModel DATA{${TestData_DIR}/LLSModel.hdf5}
  --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
  --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.fcsv
  )
set_tests_properties(BRAINSConstellationDetectorTestWarning2 PROPERTIES WILL_FAIL ON)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME LandmarksCompareTest
  COMMAND ${LAUNCH_EXE} ${LandmarksCompareEXE}
  --help
  )
endif(0)


#To define new tests for other landmark files (such as .wts files), we just need to change:
#                      TestName
#                      TEST_PROGRAM_ARG (which is the BCD flag for tested output argument)
#                      TEST_BASELINE (just add proper file to the TestData)
#                      TEST_TEMP_OUTPUT_LANDMARKS (notice to its format)
#you can use MORE_ARG and MORE_ARG_CONTENT to add new arguments to the BCD if you need to pass them for the test

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputLandmarksInACPCAlignedSpace)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${CMAKE_COMMAND}
        -D TEST_PROGRAM=${BRAINSConstellationDetectorEXE}
        -D TEST_COMPARE_PROGRAM=${LandmarksCompareEXE}
        -D TEST_BASELINE=DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv}
        -D TEST_INPUT_VOLUME=DATA{${TestData_DIR}/T1.nii.gz}
        -D TEST_INPUT_TEMPLATEMODEL=DATA{${TestData_DIR}/T1.mdl}
        -D TEST_INPUT_LLSMODEL=DATA{${TestData_DIR}/LLSModel.hdf5}
        -D TEST_TEMP_OUTPUT_RESAMPLEDVOLUME=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
        -D TEST_TEMP_OUTPUT_LANDMARKS=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.fcsv
        -D TEST_PROGRAM_ARG=outputLandmarksInACPCAlignedSpace
        -D MORE_ARG=
        -D MORE_ARG_CONTENT=
        -P ${CMAKE_CURRENT_SOURCE_DIR}/LandmarksCompareTest.cmake
        )

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputLandmarksInInputSpace)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${CMAKE_COMMAND}
         -D TEST_PROGRAM=${BRAINSConstellationDetectorEXE}
         -D TEST_COMPARE_PROGRAM=${LandmarksCompareEXE}
         -D TEST_BASELINE=DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv}
         -D TEST_INPUT_VOLUME=DATA{${TestData_DIR}/T1.nii.gz}
         -D TEST_INPUT_TEMPLATEMODEL=DATA{${TestData_DIR}/T1.mdl}
         -D TEST_INPUT_LLSMODEL=DATA{${TestData_DIR}/LLSModel.hdf5}
         -D TEST_TEMP_OUTPUT_RESAMPLEDVOLUME=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         -D TEST_TEMP_OUTPUT_LANDMARKS=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.fcsv
         -D TEST_PROGRAM_ARG=outputLandmarksInInputSpace
         -D MORE_ARG=
         -D MORE_ARG_CONTENT=
         -P ${CMAKE_CURRENT_SOURCE_DIR}/LandmarksCompareTest.cmake
         )

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputLandmarkWeights)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${CMAKE_COMMAND}
         -D TEST_PROGRAM=${BRAINSConstellationDetectorEXE}
         -D TEST_COMPARE_PROGRAM=${LandmarksCompareEXE}
         -D TEST_BASELINE=DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.wts}
         -D TEST_INPUT_VOLUME=DATA{${TestData_DIR}/T1.nii.gz}
         -D TEST_INPUT_TEMPLATEMODEL=DATA{${TestData_DIR}/T1.mdl}
         -D TEST_INPUT_LLSMODEL=DATA{${TestData_DIR}/LLSModel.hdf5}
         -D TEST_TEMP_OUTPUT_RESAMPLEDVOLUME=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         -D TEST_TEMP_OUTPUT_LANDMARKS=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.wts
         -D TEST_PROGRAM_ARG=outputLandmarkWeights
         -D MORE_ARG=
         -D MORE_ARG_CONTENT=
         -P ${CMAKE_CURRENT_SOURCE_DIR}/LandmarksCompareTest.cmake
         )

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_LandmarksPaired)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${CMAKE_COMMAND}
         -D TEST_PROGRAM=${BRAINSConstellationDetectorEXE}
         -D TEST_COMPARE_PROGRAM=${LandmarksCompareEXE}
         -D TEST_BASELINE= DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv}
         -D TEST_INPUT_VOLUME= DATA{${TestData_DIR}/T1.nii.gz}
         -D TEST_INPUT_TEMPLATEMODEL= DATA{${TestData_DIR}/T1.mdl}
         -D TEST_INPUT_LLSMODEL= DATA{${TestData_DIR}/LLSModel.hdf5}
         -D TEST_TEMP_OUTPUT_RESAMPLEDVOLUME=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         -D TEST_TEMP_OUTPUT_LANDMARKS=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.fcsv
         -D TEST_PROGRAM_ARG=outputLandmarksPaired
         -D MORE_ARG=--inputLandmarksPaired
         -D MORE_ARG_CONTENT= DATA{${TestData_DIR}/BRAINSConstellationDetectorTest_outputLandmarksInInputSpace_standard.fcsv
         -P ${CMAKE_CURRENT_SOURCE_DIR}/LandmarksCompareTest.cmake
         )

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputVerificationScript)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${CMAKE_COMMAND}
         -D TEST_PROGRAM=${BRAINSConstellationDetectorEXE}
         -D TEST_COMPARE_PROGRAM=${LandmarksCompareEXE}
         -D TEST_BASELINE= DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.sh}
         -D TEST_INPUT_VOLUME= DATA{${TestData_DIR}/T1.nii.gz}
         -D TEST_INPUT_TEMPLATEMODEL= DATA{${TestData_DIR}/T1.mdl}
         -D TEST_INPUT_LLSMODEL= DATA{${TestData_DIR}/LLSModel.hdf5}
         -D TEST_TEMP_OUTPUT_RESAMPLEDVOLUME=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         -D TEST_TEMP_OUTPUT_LANDMARKS=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.sh
         -D TEST_PROGRAM_ARG=outputVerificationScript
         -D MORE_ARG=--outputLandmarksInACPCAlignedSpace
         -D MORE_ARG_CONTENT=${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.fcsv
         -P ${CMAKE_CURRENT_SOURCE_DIR}/LandmarksCompareTest.cmake
         )

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
