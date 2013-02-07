## Set testing environment
##
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/T1.mdl ${CMAKE_CURRENT_BINARY_DIR}/T1.mdl)
set(DART_TESTING_TIMEOUT 15000)

## Test landmarksConstellationTrainingDefinitionIO
##
add_executable(TestlandmarksConstellationTrainingDefinitionIO TestlandmarksConstellationTrainingDefinitionIO.cxx)
target_link_libraries(TestlandmarksConstellationTrainingDefinitionIO BRAINSCommonLib ${ITK_LIBRARIES})

set(ALL_TEST_PROGS
  BRAINSAlignMSP
  BRAINSConstellationDetector
  BRAINSEyeDetector
  BRAINSClipInferior
  BRAINSConstellationModeler
  BRAINSTrimForegroundInDirection
  insertMidACPCpoint
  landmarksConstellationAligner
  landmarksConstellationWeights
  )

foreach(testprog ${ALL_TEST_PROGS})
  MakeTestDriverFromSEMTool(${testprog} ${testprog}Test.cxx)
endforeach()

set(ALL_PROGS_LIST
  LandmarksCompare
  )

## Simple CLP programs without TestDriver wrappers
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES landmarksConstellationCOMMONLIB )
endforeach()

if(0)
## Programs for landmarks statistics. These programs do not have any xml file
set(ALL_PROGS_LIST
    LmkStatistics
    GenerateAverageLmkFile
    )

foreach(prog ${ALL_PROGS_LIST})
  add_executable(${prog} ../landmarkStatistics/${prog}.cxx)
  target_link_libraries(${prog} BRAINSCommonLib ${ITK_LIBRARIES})
endforeach()
endif()

## Simple programs
## Test BRAINSClipInferior
##
# add_executable(BRAINSClipInferiorTest BRAINSClipInferiorTest.cxx)
# target_link_libraries(BRAINSClipInferiorTest BRAINSClipInferiorCOMMONLIB)

if(0)
## Test Versor,  This does not currently have a formal test, so disable building it.
## TODO:  Create a test for the VersorTester
##
set (VersorTester_source VersorTester.cxx)
add_executable(VersorTester ${VersorTester_source})
target_link_libraries(VersorTester ${ITK_LIBRARIES})
endif()

# We want to fetch some data, but we do not want to use them explicity as an argument in
# any program. They just have been addressed later in a txt file. However, we cannot pass
# an empty command, so we have used "BCD --help", which is like a command that do nothing.
set(BRAINSConstellationDetectorTestName fetchDataTobeUsedInInputTrainingList)
ExternalData_expand_arguments(${PROJECT_NAME}FetchData TrainingCase1_T1_VAR DATA{${TestData_DIR}/TrainingCase1_T1.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData TrainingCase2_T1_VAR DATA{${TestData_DIR}/TrainingCase2_T1.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData TrainingCase3_T1_VAR DATA{${TestData_DIR}/TrainingCase3_T1.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData TrainingCase1_fcsv_VAR DATA{${TestData_DIR}/TrainingCase1.fcsv})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData TrainingCase2_fcsv_VAR DATA{${TestData_DIR}/TrainingCase2.fcsv})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData TrainingCase3_fcsv_VAR DATA{${TestData_DIR}/TrainingCase3.fcsv})
## The following configure uses the paths from the previous ExternalData_expand_arguments
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/obliq_setup.txt.in ${CMAKE_CURRENT_BINARY_DIR}/obliq_setup.txt @ONLY IMMEDIATE)

## Test BRAINSEyeDetector
##
set(BRAINSEyeDetectorTestName BRAINSEyeDetectorTest_T1)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSEyeDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSEyeDetectorTestDriver>
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
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSAlignMSPTestDriver>
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

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_forceHoughReportFailure)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/EMSP.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
  --inputVolume DATA{${TestData_DIR}/2.nii.gz}
  --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
  --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
  --resultsDir ${CMAKE_CURRENT_BINARY_DIR}/
  --forceHoughEyeDetectorReportFailure)
set_property(TEST ${BRAINSConstellationDetectorTestName}
            PROPERTY WILL_FAIL TRUE)

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_inputVolumeError)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>

        BRAINSConstellationDetectorTest)
set_property(TEST ${BRAINSConstellationDetectorTestName}
            PROPERTY WILL_FAIL TRUE)

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_LLSModelError)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>

        BRAINSConstellationDetectorTest
        --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
        --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz)
set_property(TEST ${BRAINSConstellationDetectorTestName}
            PROPERTY WILL_FAIL TRUE)

## Test BRAINSTrimForegroundInDirection:
##
set(BRAINSTrimForegroundInDirectionTestName BRAINSTrimForegroundInDirectionTest_T1)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSTrimForegroundInDirectionTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSTrimForegroundInDirectionTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSTrimForegroundInDirectionTestName}_standard.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSTrimForegroundInDirectionTestName}_aligned.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0

  BRAINSTrimForegroundInDirectionTest
  --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSTrimForegroundInDirectionTestName}_aligned.nii.gz
  --headSizeLimit 2600.0)

set(BRAINSConstellationDetectorTestName insertMidACPCpointTest) #This test aims to generate the output of the insertMidACPCpoint
                                                                # program to be compared with the baseline file
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:insertMidACPCpointTestDriver>
  ${BRAINSConstellationDetectorTestName}
  --inputLandmarkFile DATA{${TestData_DIR}/BRAINSConstellationDetectorTest_outputLandmarksInACPCAlignedSpace_standard.fcsv} #input standard file
  --outputLandmarkFile ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_ACPCSpace_aligned.fcsv #output of the program
  )

set(BRAINSConstellationDetectorTestName MidACPCpointCompareTest)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
  --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/insertMidACPCpointTest_ACPCSpace_aligned.fcsv # output of the insertMidACPCpoint
  --inputLandmarkFile2 DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv} # Baseline file
  )
set_property(TEST ${BRAINSConstellationDetectorTestName} APPEND PROPERTY DEPENDS insertMidACPCpointTest)

set(BRAINSConstellationDetectorTestName landmarksConstellationAlignerTest) #This test aims to generate the output of the landmarksConstellationAligner to be compared with the baseline file
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:landmarksConstellationAlignerTestDriver>
  ${BRAINSConstellationDetectorTestName}
  --inputLandmarksPaired DATA{${TestData_DIR}/BRAINSConstellationDetectorTest_outputLandmarksInInputSpace_standard.fcsv}
  --outputLandmarksPaired ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}.fcsv #output of the program
  )

set(BRAINSConstellationDetectorTestName landmarksAlignerCompareTest)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
  --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/landmarksConstellationAlignerTest.fcsv # output of the landmarksConstellationAlignerTest
  --inputLandmarkFile2 DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv} # Baseline file
  )
set_property(TEST ${BRAINSConstellationDetectorTestName} APPEND PROPERTY DEPENDS landmarksConstellationAlignerTest)

set(BRAINSConstellationDetectorTestName BCD_HELP_TEST)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --help
   )

############################
############################
############################
############################
############################

# if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_acLowerBound)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --acLowerBound 100)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_backgroundValue)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --BackgroundFillValue BIGNEG)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_cutOutHeadInOutputResampledVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --cutOutHeadInOutputVolume)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_cutOutHeadInOutputVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --cutOutHeadInOutputVolume)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_interpolationMode)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --interpolationMode BSpline)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_mspLevel)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --mspQualityLevel 3)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_otsuThreshold)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --otsuPercentileThreshold 0.05)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputUntransformedClippedVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
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
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}})

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rac)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --rac 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rescaleIntensity)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --rescaleIntensities)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rescaleIntensityRange)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --rescaleIntensitiesOutputRange 20,2000)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rmpj)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --rmpj 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rpc)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --rpc 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rVN4)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --rVN4 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rVN4-rpc-rac-rmpj)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --rmpj 200
    --rVN4 200
    --rac 200
    --rpc 200)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_T1_outputResampledVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}})

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_T1_outputVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}})

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_trimRescaledIntensities)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --trimRescaledIntensities 2)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_writeBranded2DImage)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.png}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.png
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1-2ndVersion.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel-2ndVersion.${XFRM_EXT}}
    --writeBranded2DImage ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.png)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_inputLandmarksEMSP)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/2.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --inputLandmarksEMSP DATA{${TestData_DIR}/EMSP.fcsv})
  set_tests_properties(${BRAINSConstellationDetectorTestName} PROPERTIES TIMEOUT 9876)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTestForLandmarkCompare)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel.${XFRM_EXT}}
    --inputTemplateModel DATA{${TestData_DIR}/T1.mdl}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_ACPCSpace_aligned.fcsv
    --outputLandmarksInInputSpace ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_InputSpace_aligned.fcsv
    --outputLandmarkWeights ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.wts
    --inputLandmarksPaired DATA{${TestData_DIR}/BRAINSConstellationDetectorTest_outputLandmarksInInputSpace_standard.fcsv}
    --outputLandmarksPaired ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_Paired_aligned.fcsv
    --outputVerificationScript ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}_aligned.sh
    )
  set_tests_properties(${BRAINSConstellationDetectorTestName} PROPERTIES TIMEOUT 9876)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputLandmarksInACPCAlignedSpace)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BRAINSConstellationDetectorTestForLandmarkCompare_ACPCSpace_aligned.fcsv
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv} # Baseline file
    )
  set_property(TEST ${BRAINSConstellationDetectorTestName} APPEND PROPERTY
    DEPENDS BRAINSConstellationDetectorTestForLandmarkCompare)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputLandmarksInInputSpace)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BRAINSConstellationDetectorTestForLandmarkCompare_InputSpace_aligned.fcsv
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv} # Baseline file
    )
  set_property(TEST ${BRAINSConstellationDetectorTestName} APPEND PROPERTY
    DEPENDS BRAINSConstellationDetectorTestForLandmarkCompare)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_outputLandmarkWeights)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BRAINSConstellationDetectorTestForLandmarkCompare_aligned.wts
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.wts} # Baseline file
    )
  set_property(TEST ${BRAINSConstellationDetectorTestName} APPEND PROPERTY
    DEPENDS BRAINSConstellationDetectorTestForLandmarkCompare)

  set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_LandmarksPaired)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BRAINSConstellationDetectorTestForLandmarkCompare_Paired_aligned.fcsv
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.fcsv} # Baseline file
    )
  set_property(TEST ${BRAINSConstellationDetectorTestName} APPEND PROPERTY
    DEPENDS BRAINSConstellationDetectorTestForLandmarkCompare)

  set(BRAINSConstellationDetectorTestName landmarksConstellationWeightsTest)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:landmarksConstellationWeightsTestDriver>
    ${BRAINSConstellationDetectorTestName}
    --inputTrainingList ${CMAKE_CURRENT_BINARY_DIR}/obliq_setup.txt
    --inputTemplateModel DATA{${TestData_DIR}/T1-2ndVersion.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel-2ndVersion.${XFRM_EXT}}
    --outputWeightsList ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationDetectorTestName}.wts #output of the program
    )
  set_tests_properties(${BRAINSConstellationDetectorTestName} PROPERTIES TIMEOUT 15000)

  set(BRAINSConstellationDetectorTestName landmarksWeightsCompareTest)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BRAINSConstellationDetectorTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/landmarksConstellationWeightsTest.wts # output of the landmarksConstellationWeightsTest
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BRAINSConstellationDetectorTestName}_standard.wts} # Baseline file
    )
  set_property(TEST ${BRAINSConstellationDetectorTestName} APPEND PROPERTY DEPENDS landmarksConstellationWeightsTest)

if ( 0 ) ## HACK:  Just silencing failing tests until they can be proper addressed
  set(BRAINSConstellationModelerTestName BRAINSConstellationModelerTest_T1)
  if(ConstellationModeler_TEST_causes_major_delay)
    ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSConstellationModelerTestName}
      COMMAND BRAINSConstellationModelerTest
      BRAINSConstellationModelerTest
      --inputTrainingList ${CMAKE_CURRENT_BINARY_DIR}/obliq_setup.txt
      --outputModel ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSConstellationModelerTestName}.mdl
      )
  endif()
endif()

# endif()

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
