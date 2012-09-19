
MakeTestDriverFromSEMTool(BRAINSROIAuto BRAINSROIAutoTest.cxx)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSROIAutoTest_GenerateBrain_seg
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSROIAutoTestDriver>
  --compare DATA{${TestData_DIR}/BRAINSROIAutoTest_GenerateBrainMask.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSROIAutoTest_GenerateBrain_seg.nii.gz
  BRAINSROIAutoTest
  --inputVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputROIMaskVolume ${CMAKE_CURRENT_BINARY_DIR}/BRAINSROIAutoTest_GenerateBrain_seg.nii.gz
  --ROIAutoDilateSize 10
  --outputVolumePixelType short
  --maskOutput
  --cropOutput
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/BRAINSROIAutoTest_GenerateBrain_clipped_cropped.nii.gz
  )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
