
MakeTestDriverFromSEMTool(BRAINSMush BRAINSMushTest.cxx)

configure_file(${BRAINSMush_SOURCE_DIR}/CTestCustom.ctest ${BRAINSMush_BINARY_DIR}/CTestCustom.ctest COPYONLY)

##  Test goes here:
set(BRAINSMushTestName BRAINSMushTest_MaskGiven)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSMushTestName} COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSMushTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSMushTestName}_MushResult_BrainMask.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSMushTest
  --inputFirstVolume DATA{${TestData_DIR}/T1Aligned.nii.gz}
  --inputSecondVolume DATA{${TestData_DIR}/T2Aligned.nii.gz}
  --inputMaskVolume DATA{${TestData_DIR}/BrainMaskInit.nii.gz}
  --outputWeightsFile ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_weights.txt
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_MushTest.nii.gz
  --outputMask ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
)

set(BRAINSMushTestName BRAINSMushTest_NoMaskGiven)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSMushTestName} COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSMushTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSMushTestName}_MushResult_BrainMask.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSMushTest
  --inputFirstVolume DATA{${TestData_DIR}/T1Aligned.nii.gz}
  --inputSecondVolume DATA{${TestData_DIR}/T2Aligned.nii.gz}
  --outputWeightsFile ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_weights.txt
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_MushTest.nii.gz
  --outputMask ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSMushTestName}_MushTest_BrainMask.nii.gz
)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
