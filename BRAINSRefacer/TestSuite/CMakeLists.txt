#TODO: ADD overview of tests


#TODO: create test data that can be shared
set(RefacerTestData_Input_DIR "$ENV{HOME}/testDataNEW/BRAINSRefacerData")
set(RefacerTestData_Output_DIR "$ENV{HOME}/testDataNEW/Output")
message(WARNING "TestData should be placed into directory ${RefacerTestData_DIR}")

### setup test programs/drivers

MakeTestDriverFromSEMTool(BRAINSRefacer BRAINSRefacerTest.cxx)

#add all simple CLP programs to this list
set(ALL_PROGS_LIST_BRAINSDEFACETEST
  BRAINSRefacerBrainDataTest
  BRAINSRefacerCheckDeformationAreaLabelMap
  BRAINSRefacerCheckDeformationAreaLandmarks
  )

#For simple CLP programs without TestDriver wrappers
foreach(prog ${ALL_PROGS_LIST_BRAINSDEFACETEST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

### setup pixel types to be tested on

# component types from itkImageIObase - 2016/7/8
# UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE
set(PIXEL_LIST_BRAINSREFACER
  UCHAR
  CHAR
  USHORT
  SHORT
  UINT
  INT
  ULONG
  LONG
  FLOAT
  DOUBLE
  )

### Run all the tests against each pixel type

# Loop through all pixel types for each test
foreach(PTYPE ${PIXEL_LIST_BRAINSREFACER} )
  set(BRR_inPrefix "${RefacerTestData_Input_DIR}/${PTYPE}/${PTYPE}_")
  set(BRR_outPrefix "${RefacerTestData_Output_DIR}/${PTYPE}_")

  # basic test to make sure landmarks based deformation will run
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_${PTYPE}_RunTest_Landmarks
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    BRAINSRefacerTest
    --inputImage
    ${BRR_inPrefix}Subject.nii.gz
    --landmarks
    ${BRR_inPrefix}Landmarks.fcsv
    --deformedImageName
    ${BRR_outPrefix}LandmarksRefaced_Subject.nii.gz
    )

    # check that no brain data has changed in the landmarks based deformation
    # note this test requires a reference/ground truth brain label map
    set(PARENT_TEST BRAINSRefacerTest_${PTYPE}_RunTest_Landmarks)
    set(RefacerTestName BRAINSRefacerTest_${PTYPE}_CheckBrainData_Landmarks)
    ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
      COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
      --inputOriginal
      ${BRR_inPrefix}Subject.nii.gz
      --inputRefaced
      ${BRR_outPrefix}LandmarksRefaced_Subject.nii.gz
      --labelmap
      ${BRR_inPrefix}BrainLabels.nii.gz
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  # check to see that non brain data changed for landmarks
  set(PARENT_TEST BRAINSRefacerTest_${PTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerTest_${PTYPE}_CheckDeformationArea_Landmarks)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLandmarks>
    --inputOriginal
    ${BRR_inPrefix}Subject.nii.gz
    --inputRefaced
    ${BRR_outPrefix}LandmarksRefaced_Subject.nii.gz
    --landmarks
    ${BRR_inPrefix}Landmarks.fcsv
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  set(PARENT_TEST BRAINSRefacerTest_${PTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerTest_${PTYPE}_CheckNonDeformationArea_Landmarks)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLandmarks>
    --inputOriginal
    ${BRR_inPrefix}Subject.nii.gz
    --inputRefaced
    ${BRR_outPrefix}LandmarksRefaced_Subject.nii.gz
    --landmarks
    ${BRR_inPrefix}Landmarks.fcsv
    --checkNonDeformedArea
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})


  ####

  # basic test to make sure brain labelmap based deformation will run
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_${PTYPE}_RunTest_Labelmap
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    BRAINSRefacerTest
    --inputImage
    ${BRR_inPrefix}Subject.nii.gz
    --labelmap
    ${BRR_inPrefix}BrainLabels.nii
    --deformedImageName
    ${BRR_outPrefix}LabelMap_RefacedSubject.nii.gz
    --labelmapMask
    )

  #check that no brain data has been changed in labelmap based deformation
  # note this test requires a reference/ground truth brain label map
  set(PARENT_TEST BRAINSRefacerTest_${PTYPE}_RunTest_Labelmap)
  set(RefacerTestName BRAINSRefacerTest_${PTYPE}_CheckBrainData_Labelmap)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
    --inputOriginal
    ${BRR_inPrefix}Subject.nii.gz
    --inputRefaced
    ${BRR_outPrefix}LabelMap_RefacedSubject.nii.gz
    --labelmap
    ${BRR_inPrefix}BrainLabels.nii.gz
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})


  # check to see that non brain data has changed for labelmap
  set(PARENT_TEST BRAINSRefacerTest_${PTYPE}_RunTest_Labelmap)
  set(RefacerTestName BRAINSRefacerTest_${PTYPE}_CheckDeformationArea_Labelmap)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLabelMap>
    --inputOriginal
    ${BRR_inPrefix}Subject.nii.gz
    --inputRefaced
    ${BRR_outPrefix}LabelMap_RefacedSubject.nii.gz
    --labelmap
    ${BRR_inPrefix}BrainLabels.nii
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  set(PARENT_TEST BRAINSRefacerTest_${PTYPE}_RunTest_Labelmap)
  set(RefacerTestName BRAINSRefacerTest_${PTYPE}_CheckNonDeformationArea_Labelmap)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLabelMap>
    --inputOriginal
    ${BRR_inPrefix}Subject.nii.gz
    --inputRefaced
    ${BRR_outPrefix}LabelMap_RefacedSubject.nii.gz
    --labelmap
    ${BRR_inPrefix}BrainLabels.nii
    --checkNonDeformedArea
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #TODO: Test to reload bSpline and apply it to the same image, and compare to see if the output is the same

  #TODO: Test to reload bSpline and apply it to a new image

  #TODO: Test that the output pixel type is the same as the input

  #TODO: Test that the deformation field is symmetrical

endforeach()

###

#TODO: once the above TODOs are all finished and checked, remove the tests below
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
