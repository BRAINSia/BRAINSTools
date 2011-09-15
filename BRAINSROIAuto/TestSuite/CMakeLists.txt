
MakeTestDriverFromSEMTool(BRAINSROIAuto BRAINSROIAutoTest.cxx)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSROIAutoTest_GenerateBrainMask
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSROIAutoTestDriver>
  --compare DATA{${TestData_DIR}/BRAINSROIAutoTest_GenerateBrainMask.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSROIAutoTest_GenerateBrainMask.nii.gz
  BRAINSROIAutoTest
  --inputVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputROIMaskVolume ${CMAKE_CURRENT_BINARY_DIR}/BRAINSROIAutoTest_GenerateBrainMask.nii.gz
  --ROIAutoDilateSize 10
  --outputVolumePixelType uchar
  )
