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
  --labelmap
    ${RefacerTestData_DIR}/JointFusion_HDAtlas20_2015_dustCleaned_label.nii.gz
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

set(PARENT_TEST BRAINSRefacerTest_RunTest)
set(RefacerTestName BRAINSRefacerBrainDataTest})
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
        COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
        )
set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})