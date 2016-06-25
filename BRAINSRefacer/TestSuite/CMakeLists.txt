MakeTestDriverFromSEMTool(BRAINSRefacer BRAINSRefacerTest.cxx)
set(RefacerTestData_DIR "/scratch/testData/BRAINSRefacerData")
message(WARNING "TestData should be placed into directory ${RefacerTestData_DIR}")


## add all sipmle CLP programs to cmakelist
set(ALL_PROGS_LIST_BRAINSDEFACETEST
        BRAINSRefacerBrainDataTest
        )

## Simple CLP programs without TestDriver wrappers
foreach(prog ${ALL_PROGS_LIST_BRAINSDEFACETEST})
    StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_RunTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
  BRAINSRefacerTest
  --inputImage
    ${RefacerTestData_DIR}/t1_average_BRAINSABC.nii.gz
  --landmarks
    ${RefacerTestData_DIR}/BCD_ACPC_Landmarks.fcsv
  --deformedImageName
    ${RefacerTestData_DIR}/defacedSubject.nii.gz
  --outputMask
    ${RefacerTestData_DIR}/outputMask.nii.gz
  --distanceMapFileName
    ${RefacerTestData_DIR}/distanceMap.nii.gz
  --bSplineFileName
    ${RefacerTestData_DIR}/bSpline.h5
  --smoothDisplacementName
    ${RefacerTestData_DIR}/bSplineTimesDistanceMap.nii.gz
  --diffImageName
    ${RefacerTestData_DIR}/diffImage.nii.gz
  --finalTransformFileName
    ${RefacerTestData_DIR}/finalTransform.h5
  )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_RunTest_Labelmap
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
  BRAINSRefacerTest
  --inputImage
    ${RefacerTestData_DIR}/t1_average_BRAINSABC.nii.gz
  --labelmap
    ${RefacerTestData_DIR}/JointFusion_HDAtlas20_2015_dustCleaned_label.nii
  --deformedImageName
    ${RefacerTestData_DIR}/LabelMap_defacedSubject.nii.gz
  --outputMask
    ${RefacerTestData_DIR}/LabelMap_outputMask.nii.gz
  --distanceMapFileName
    ${RefacerTestData_DIR}/LabelMap_distanceMap.nii.gz
  --bSplineFileName
    ${RefacerTestData_DIR}/LabelMap_bSpline.h5
  --smoothDisplacementName
    ${RefacerTestData_DIR}/LabelMap_bSplineTimesDistanceMap.nii.gz
  --diffImageName
    ${RefacerTestData_DIR}/LabelMap_diffImage.nii.gz
  --finalTransformFileName
    ${RefacerTestData_DIR}/LabelMap_finalTransform.h5
  --labelmapMask
  )
set(PARENT_TEST BRAINSRefacerTest_RunTest_Labelmap)
set(RefacerTestName BRAINSRefacerBrainDataTest_LabelMapMask)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
        COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
        --inputOriginal ${RefacerTestData_DIR}/t1_average_BRAINSABC.nii.gz
        --inputRefaced ${RefacerTestData_DIR}/LabelMap_defacedSubject.nii.gz
        --labelmap ${RefacerTestData_DIR}/JointFusion_HDAtlas20_2015_dustCleaned_label.nii.gz
        --labelmapMask
  )
set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

set(PARENT_TEST BRAINSRefacerTest_RunTest)
set(RefacerTestName BRAINSRefacerBrainDataTest_LandmarkMapMask)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
  --inputOriginal ${RefacerTestData_DIR}/t1_average_BRAINSABC.nii.gz
  --inputRefaced ${RefacerTestData_DIR}/defacedSubject.nii.gz
  --landmarks ${RefacerTestData_DIR}/BCD_ACPC_Landmarks.fcsv
  )
set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

# check to see that non brain data changed
set(PARENT_TEST BRAINSRefacerTest_RunTest)
set(RefacerTestName BRAINSRefacerBrainDataTest_didFaceChange)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
        COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
        --inputOriginal ${RefacerTestData_DIR}/t1_average_BRAINSABC.nii.gz
        --inputRefaced ${RefacerTestData_DIR}/defacedSubject.nii.gz
        --landmarks ${RefacerTestData_DIR}/BCD_ACPC_Landmarks.fcsv
        --checkNonBrainData
        )
set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})


#try running with different sets of data
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_RunTest2
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
  BRAINSRefacerTest
  --inputImage
  ${RefacerTestData_DIR}/DONOTSHARE/0005/t1_average_BRAINSABC.nii.gz
  --labelmap
  ${RefacerTestData_DIR}/DONOTSHARE/0005/complete_brainlabels_seg.nii.gz
  --landmarks

  ${RefacerTestData_DIR}/DONOTSHARE/0005/BCD_ACPC_Landmarks.fcsv
  --deformedImageName
  ${RefacerTestData_DIR}/RunTest2_defacedSubject.nii.gz
  --outputMask
  ${RefacerTestData_DIR}/RunTest2_outputMask.nii.gz
  --distanceMapFileName
  ${RefacerTestData_DIR}/RunTest2_distanceMap.nii.gz
  --bSplineFileName
  ${RefacerTestData_DIR}/RunTest2_bSpline.h5
  --smoothDisplacementName
  ${RefacerTestData_DIR}/RunTest2_bSplineTimesDistanceMap.nii.gz
  --diffImageName
  ${RefacerTestData_DIR}/RunTest2_diffImage.nii.gz
  --finalTransformFileName
  ${RefacerTestData_DIR}/RunTest2_finalTransform.h5
  --labelmapMask
  )


ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_RunTest3_labelmap
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
  BRAINSRefacerTest
  --inputImage
  ${RefacerTestData_DIR}/DONOTSHARE/0006/t2_average_BRAINSABC.nii.gz
  --landmarks
  ${RefacerTestData_DIR}/DONOTSHARE/0006/BCD_ACPC_Landmarks.fcsv
  --deformedImageName
  ${RefacerTestData_DIR}/RunTest3_defacedSubject.nii.gz
  --outputMask
  ${RefacerTestData_DIR}/RunTest3_outputMask.nii.gz
  --distanceMapFileName
  ${RefacerTestData_DIR}/RunTest3_distanceMap.nii.gz
  --bSplineFileName
  ${RefacerTestData_DIR}/RunTest3_bSpline.h5
  --smoothDisplacementName
  ${RefacerTestData_DIR}/RunTest3_bSplineTimesDistanceMap.nii.gz
  --diffImageName
  ${RefacerTestData_DIR}/RunTest3_diffImage.nii.gz
  --finalTransformFileName
  ${RefacerTestData_DIR}/RunTest3_finalTransform.h5
  --labelmapMask
  --labelmap
  ${RefacerTestData_DIR}/DONOTSHARE/0006/complete_brainlabels_seg.nii.gz
  )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_RunTest3_landmark
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
  BRAINSRefacerTest
  --inputImage
  ${RefacerTestData_DIR}/DONOTSHARE/0006/t2_average_BRAINSABC.nii.gz
  --landmarks
  ${RefacerTestData_DIR}/DONOTSHARE/0006/BCD_ACPC_Landmarks.fcsv
  --deformedImageName
  ${RefacerTestData_DIR}/RunTest3_landmark_defacedSubject.nii.gz
  --outputMask
  ${RefacerTestData_DIR}/RunTest3_landmark_outputMask.nii.gz
  --distanceMapFileName
  ${RefacerTestData_DIR}/RunTest3_landmark_distanceMap.nii.gz
  --bSplineFileName
  ${RefacerTestData_DIR}/RunTest3_landmark_bSpline.h5
  --smoothDisplacementName
  ${RefacerTestData_DIR}/RunTest3_landmark_bSplineTimesDistanceMap.nii.gz
  --diffImageName
  ${RefacerTestData_DIR}/RunTest3_landmark_diffImage.nii.gz
  --finalTransformFileName
  ${RefacerTestData_DIR}/RunTest3_landmark_finalTransform.h5

  )

#try running with different sets of data
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_RunTest4
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
  BRAINSRefacerTest
  --inputImage
  ${RefacerTestData_DIR}/DONOTSHARE/0007/t1_average_BRAINSABC.nii.gz
  --landmarks
  ${RefacerTestData_DIR}/DONOTSHARE/0007/BCD_ACPC_Landmarks.fcsv
  --deformedImageName
  ${RefacerTestData_DIR}/RunTest4_defacedSubject.nii.gz
  --outputMask
  ${RefacerTestData_DIR}/RunTest4_outputMask.nii.gz
  --distanceMapFileName
  ${RefacerTestData_DIR}/RunTest4_distanceMap.nii.gz
  --bSplineFileName
  ${RefacerTestData_DIR}/RunTest4_bSpline.h5
  --smoothDisplacementName
  ${RefacerTestData_DIR}/RunTest4_bSplineTimesDistanceMap.nii.gz
  --diffImageName
  ${RefacerTestData_DIR}/RunTest4_diffImage.nii.gz
  --finalTransformFileName
  ${RefacerTestData_DIR}/RunTest4_finalTransform.h5
  --labelmapMask
  --labelmap
  ${RefacerTestData_DIR}/DONOTSHARE/0007/complete_brainlabels_seg.nii.gz
  )



#try running with different sets of data
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_RunTest5
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
  BRAINSRefacerTest
  --inputImage
  ${RefacerTestData_DIR}/DONOTSHARE/0007/t2_average_BRAINSABC.nii.gz
  --landmarks
  ${RefacerTestData_DIR}/DONOTSHARE/0007/BCD_ACPC_Landmarks.fcsv
  --deformedImageName
  ${RefacerTestData_DIR}/RunTest5_defacedSubject.nii.gz
  --outputMask
  ${RefacerTestData_DIR}/RunTest5_outputMask.nii.gz
  --distanceMapFileName
  ${RefacerTestData_DIR}/RunTest5_distanceMap.nii.gz
  --bSplineFileName
  ${RefacerTestData_DIR}/RunTest5_bSpline.h5
  --smoothDisplacementName
  ${RefacerTestData_DIR}/RunTest5_bSplineTimesDistanceMap.nii.gz
  --diffImageName
  ${RefacerTestData_DIR}/RunTest5_diffImage.nii.gz
  --finalTransformFileName
  ${RefacerTestData_DIR}/RunTest5_finalTransform.h5
  )
