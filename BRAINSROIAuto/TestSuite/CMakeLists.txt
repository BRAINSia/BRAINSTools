
MakeTestDriverFromSEMTool(BRAINSROIAuto BRAINSROIAutoTest.cxx)

if(BRAINS_BUILD)
  set(TESTING_BINARY_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
else(BRAINS_BUILD)
  set(TESTING_BINARY_DIR ${BRAINSROIAuto_BINARY_DIR}/bin)
endif(BRAINS_BUILD)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSROIAutoTest_GenerateBrainMask
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSROIAutoTestDriver>
  --compare DATA{${TestData_DIR}/BRAINSROIAutoTest_GenerateBrainMask.result.nii.gz}
                 ${TESTING_BINARY_DIR}/BRAINSROIAutoTest_GenerateBrainMask.nii.gz
  --inputVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputROIMaskVolume ${TESTING_BINARY_DIR}/BRAINSROIAutoTest_GenerateBrainMask.nii.gz
  --ROIAutoDilateSize 10
  --outputVolumePixelType uchar
  )
