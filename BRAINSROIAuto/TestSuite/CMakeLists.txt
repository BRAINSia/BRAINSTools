
MakeTestDriverFromSEMTool(BRAINSROIAuto BRAINSROIAutoTest.cxx)

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME BRAINSROIAutoTest_GenerateBrain_seg
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

# Declare per-test CPU consumption so CTest manages parallelism natively.
get_property(_all_tests DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
set_tests_properties(${_all_tests} PROPERTIES PROCESSORS 3)
