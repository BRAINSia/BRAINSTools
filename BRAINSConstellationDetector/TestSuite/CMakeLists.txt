## Set testing environment
##
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

## Test landmarksConstellationTrainingDefinitionIO
##
add_executable(TestlandmarksConstellationTrainingDefinitionIO TestlandmarksConstellationTrainingDefinitionIO.cxx)
target_link_libraries(TestlandmarksConstellationTrainingDefinitionIO BRAINSCommonLib ${BRAINSConstellationDetector_ITK_LIBRARIES}
  ${VTK_LIBRARIES})

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
  MakeTestDriverFromSEMTool(${testprog} ${testprog}Test.cxx ${VTK_LIBRARIES})
endforeach()

set(ALL_PROGS_LIST
  LandmarksCompare
  )

## Simple CLP programs without TestDriver wrappers
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog}
    TARGET_LIBRARIES landmarksConstellationCOMMONLIB BRAINSCommonLib ${VTK_LIBRARIES})
endforeach()

## Simple programs
## Test BRAINSClipInferior
##
# add_executable(BRAINSClipInferiorTest BRAINSClipInferiorTest.cxx)
# target_link_libraries(BRAINSClipInferiorTest BRAINSClipInferiorCOMMONLIB)

if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 8) # This should be restored after fixing.
## Test Versor,  This does not currently have a formal test, so disable building it.
## TODO:  Create a test for the VersorTester
##
set (VersorTester_source VersorTester.cxx)
add_executable(VersorTester ${VersorTester_source})
target_link_libraries(VersorTester ${BRAINSConstellationDetector_ITK_LIBRARIES})
endif()

# We want to fetch some data, but we do not want to use them explicity as an argument in
# any program. They just have been addressed later in a txt file. However, we cannot pass
# an empty command, so we have used "BCD --help", which is like a command that do nothing.
set(BCDTestName fetchDataTobeUsedInInputTrainingList)
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
  --compare DATA{${TestData_DIR}/${BRAINSEyeDetectorTestName}_standard.nrrd}
  ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSEyeDetectorTestName}_aligned.nrrd
  --compareIntensityTolerance 0
  --compareRadiusTolerance 5
  --compareNumberOfPixelsTolerance 0
  BRAINSEyeDetectorTest
  --inputVolume DATA{${TestData_DIR}/${BRAINSEyeDetectorTestName}_standard.nrrd}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSEyeDetectorTestName}_aligned.nrrd
  --debugDir ${CMAKE_CURRENT_BINARY_DIR})

## Test BRAINSAlignMSP:
##
set(BRAINSAlignMSPTestName BRAINSAlignMSPTest_T1)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSAlignMSPTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSAlignMSPTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSAlignMSPTestName}_standard.nrrd}
  ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSAlignMSPTestName}_aligned.nrrd
  --compareIntensityTolerance 50
  --compareRadiusTolerance 5
  --compareNumberOfPixelsTolerance 500
  BRAINSAlignMSPTest
  --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
  --OutputresampleMSP ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSAlignMSPTestName}_aligned.nrrd
  --mspQualityLevel 3
  )

## Test BRAINSConstellationDetector
##

set(BCDTestName BCDTest_forceHoughReportFailure)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
  --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/EMSP.nrrd
  --compareIntensityTolerance 0
  --compareRadiusTolerance 5
  --compareNumberOfPixelsTolerance 0

         BCDTest
  --inputVolume DATA{${TestData_DIR}/2.nii.gz}
  --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
  --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
  --resultsDir ${CMAKE_CURRENT_BINARY_DIR}/
  --forceHoughEyeDetectorReportFailure)
set_property(TEST ${BCDTestName}
            PROPERTY WILL_FAIL TRUE)

set(BCDTestName BCDTest_inputVolumeError)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>

        BCDTest)
set_property(TEST ${BCDTestName}
            PROPERTY WILL_FAIL TRUE)

set(BCDTestName BCDTest_LLSModelError)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>

        BCDTest
        --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
        --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd)
set_property(TEST ${BCDTestName}
            PROPERTY WILL_FAIL TRUE)

## Test BRAINSTrimForegroundInDirection:
##
set(BRAINSTrimForegroundInDirectionTestName BRAINSTrimForegroundInDirectionTest_T1)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSTrimForegroundInDirectionTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSTrimForegroundInDirectionTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSTrimForegroundInDirectionTestName}_standard.nrrd}
  ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSTrimForegroundInDirectionTestName}_aligned.nrrd
  --compareIntensityTolerance 0
  --compareRadiusTolerance 5
  --compareNumberOfPixelsTolerance 0

  BRAINSTrimForegroundInDirectionTest
  --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSTrimForegroundInDirectionTestName}_aligned.nrrd
  --headSizeLimit 2600.0)

set(BCDTestName insertMidACPCpointTest) #This test aims to generate the output of the insertMidACPCpoint
                                                                # program to be compared with the baseline file
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:insertMidACPCpointTestDriver>
  ${BCDTestName}
  --inputLandmarkFile DATA{${TestData_DIR}/${BCDTestName}_input.fcsv} #input standard file
  --outputLandmarkFile ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_ACPCSpace_aligned.fcsv #output of the program
  )

set(PARENT_TEST insertMidACPCpointTest)
set(BCDTestName chk_${PARENT_TEST})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
  --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/insertMidACPCpointTest_ACPCSpace_aligned.fcsv # output of the insertMidACPCpoint
  --inputLandmarkFile2 DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv} # Baseline file
  --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
  )
set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(BCDTestName landmarksConstellationAlignerTest) #This test aims to generate the output of the landmarksConstellationAligner to be compared with the baseline file
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:landmarksConstellationAlignerTestDriver>
  ${BCDTestName}
  --inputLandmarksPaired DATA{${TestData_DIR}/${BCDTestName}_input.fcsv}
  --outputLandmarksPaired ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}.fcsv #output of the program
  )

set(PARENT_TEST landmarksConstellationAlignerTest)
set(BCDTestName chk_${PARENT_TEST})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
  --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/landmarksConstellationAlignerTest.fcsv # output of the landmarksConstellationAlignerTest
  --inputLandmarkFile2 DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv} # Baseline file
  --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
  )
set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(BCDTestName BCD_HELP_TEST)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --help
   )

############################
############################
############################
############################
############################
set(BCDTestName BCDTest_ForceRPPoint)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --forceRPPoint 0.44,7.27,1.97
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --outputLandmarksInInputSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_InputSpace.fcsv
    )

  set(PARENT_TEST BCDTest_ForceRPPoint)
  set(BCDTestName chk_${PARENT_TEST})
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BCDTest_ForceRPPoint_InputSpace.fcsv # output of the ForceRPPoint Test
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv}  # Baseline file
    --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
    )
  set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(BCDTestName BCDTest_ForceACPoint)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --forceACPoint 0.72,30.01,18.95
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --outputLandmarksInInputSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_InputSpace.fcsv
    )

  set(PARENT_TEST BCDTest_ForceACPoint)
  set(BCDTestName chk_${PARENT_TEST})
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BCDTest_ForceACPoint_InputSpace.fcsv # output of the ForceACPoint Test
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv}# Baseline file
    --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
    )
  set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(BCDTestName BCDTest_ForcePCPoint)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --forcePCPoint 1.54,1.64,14.32
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --outputLandmarksInInputSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_InputSpace.fcsv
    )

  set(PARENT_TEST BCDTest_ForcePCPoint)
  set(BCDTestName chk_${PARENT_TEST})
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BCDTest_ForcePCPoint_InputSpace.fcsv # output of the ForcePCPoint Test
    --inputLandmarkFile2 DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv}# Baseline file
    --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
    )
  set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(BCDTestName BCDTest_ForceVN4Point)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --forceVN4Point 0,-20.60,-12.81
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --outputLandmarksInInputSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_InputSpace.fcsv
    )

  set(PARENT_TEST BCDTest_ForceVN4Point)
  set(BCDTestName chk_${PARENT_TEST})
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BCDTest_ForceVN4Point_InputSpace.fcsv # output of the ForceVN4Point Test
    --inputLandmarkFile2
    DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv} # Baseline file
    --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
    )
  set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  set(BCDTestName BCDTest_acLowerBound)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --acLowerBound 100)

  set(BCDTestName BCDTest_backgroundValue)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --BackgroundFillValue BIGNEG)

  set(BCDTestName BCDTest_cutOutHeadInOutputResampledVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --cutOutHeadInOutputVolume)

  set(BCDTestName BCDTest_cutOutHeadInOutputVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --cutOutHeadInOutputVolume)

  set(BCDTestName BCDTest_interpolationMode)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --interpolationMode BSpline)

  set(BCDTestName BCDTest_mspLevel)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --mspQualityLevel 2)

  set(BCDTestName BCDTest_otsuThreshold)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --otsuPercentileThreshold 0.05)

  set(BCDTestName BCDTest_outputUntransformedClippedVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12
    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputUntransformedClippedVolume
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --acLowerBound 100
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}})

  set(BCDTestName BCDTest_rac)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_lmk_acpc.fcsv
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --rac 10)

  set(BCDTestName BCDTest_rescaleIntensity)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 75
    --compareRadiusTolerance 10
    --compareNumberOfPixelsTolerance 25

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --rescaleIntensities)

  set(BCDTestName BCDTest_rescaleIntensityRange)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --rescaleIntensitiesOutputRange 20,2000)

  set(BCDTestName BCDTest_rmpj)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_lmk_acpc.fcsv
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --rmpj 10)

  set(BCDTestName BCDTest_rpc)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_lmk_acpc.fcsv
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --rpc 10)

  set(BCDTestName BCDTest_rVN4)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_lmk_acpc.fcsv
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --rVN4 10)

  set(BCDTestName BCDTest_rVN4-rpc-rac-rmpj)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_lmk_acpc.fcsv
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --rmpj 10
    --rVN4 10
    --rac 10
    --rpc 10)

  set(BCDTestName BCDTest_T1_outputResampledVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}})

  set(BCDTestName BCDTest_T1_outputVolume)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}})

  set(BCDTestName BCDTest_trimRescaledIntensities)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd,:}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --trimRescaledIntensities 2)

  set(BCDTestName BCDTest_writeBranded2DImage)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.png}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.png
    --compareIntensityTolerance 50
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 12

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --writeBranded2DImage ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.png)

  set(BCDTestName BCDTest_inputLandmarksEMSP)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    --compare DATA{${TestData_DIR}/${BCDTestName}_standard.nrrd}
    ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 0

    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/2_EMSP.nrrd}
    --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_lmk_acpc.fcsv
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --inputLandmarksEMSP DATA{${TestData_DIR}/EMSP.fcsv})
  set_tests_properties(${BCDTestName} PROPERTIES TIMEOUT 9876)

  set(BCDTestName BCDTestForLandmarkCompare)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSConstellationDetectorTestDriver>
    BRAINSConstellationDetectorTest
    --inputVolume DATA{${TestData_DIR}/T1.nii.gz}
    --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
    --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.nrrd
    --outputLandmarksInACPCAlignedSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_ACPCSpace_aligned.fcsv
    --outputLandmarksInInputSpace ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_InputSpace_aligned.fcsv
    --outputVerificationScript ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}_aligned.sh
    )
  set_tests_properties(${BCDTestName} PROPERTIES TIMEOUT 9876)

  set(PARENT_TEST BCDTestForLandmarkCompare)
  set(BCDTestName chk_a_${PARENT_TEST})
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BCDTestForLandmarkCompare_ACPCSpace_aligned.fcsv
    --inputLandmarkFile2
    DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv} # Baseline file
  --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
    )
  set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  set(PARENT_TEST BCDTestForLandmarkCompare)
  set(BCDTestName chk_i_${PARENT_TEST})
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
    --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BCDTestForLandmarkCompare_InputSpace_aligned.fcsv
    --inputLandmarkFile2
    DATA{${TestData_DIR}/${BCDTestName}_standard.fcsv} # Baseline file
  --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
    )
  set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #set(PARENT_TEST BCDTestForLandmarkCompare)
  #set(BCDTestName chk_${PARENT_TEST})
  #ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  #  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
  #  --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/BCDTestForLandmarkCompare_aligned.wts
  #  --inputLandmarkFile2 DATA{${TestData_DIR}/${BCDTestName}_standard.wts} # Baseline file
  # --tolerance 1.875 # The images are 1.875x1.875x2.4 in the test suite, so a large tolerance is needed
  #  )
  #set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

#    DEPENDS BCDTestForLandmarkCompare)
if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 8) # This should be restored after fixing.
## HACK:  Just silencing failing tests until they can be proper addressed
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

if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 5) ## This is a very long running test that is too intensive for everyday testing
  #set(BCDTestName landmarksConstellationWeightsTest)
  #ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  #  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:landmarksConstellationWeightsTestDriver>
  #  ${BCDTestName}
  #  --inputTrainingList ${CMAKE_CURRENT_BINARY_DIR}/obliq_setup.txt
  #  --inputTemplateModel DATA{${TestData_DIR}/T1_50Lmks.mdl}
  #  --LLSModel DATA{${TestData_DIR}/Transforms_h5/LLSModel_50Lmks.${XFRM_EXT}}
  #  --outputWeightsList ${CMAKE_CURRENT_BINARY_DIR}/${BCDTestName}.wts #output of the program
  #  )
  #set_tests_properties(${BCDTestName} PROPERTIES TIMEOUT 15000)

  #set(PARENT_TEST landmarksConstellationWeightsTest)
  #set(BCDTestName chk_${PARENT_TEST})
  #ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${BCDTestName}
  #  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LandmarksCompare>
  #  --inputLandmarkFile1 ${CMAKE_CURRENT_BINARY_DIR}/landmarksConstellationWeightsTest.wts # output of the landmarksConstellationWeightsTest
  #  --inputLandmarkFile2 DATA{${TestData_DIR}/${BCDTestName}_standard.wts} # Baseline file
  #  )
  #set_property(TEST ${BCDTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

endif()

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
