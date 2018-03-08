###  This is the test suite

MakeTestDriverFromSEMTool(BRAINSContinuousClass BRAINSContinuousClassTest.cxx)
add_dependencies(BRAINSContinuousClassTestDriver InstallReferenceAtlas) ## Needed to ensure data is installed

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME ValidateBRAINSContinuousClassTest_nii
  COMMAND ${LAUNCH_EXE} ${CMAKE_BIN_DIR}/BRAINSContinuousClassTestDriver
  --compare
  DATA{${TestData_DIR}/BRAINSContinuousClassTest_baseline.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSContinuousClassTest_test.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSContinuousClassTest
  --inputT1Volume ${ATLAS_INSTALL_DIRECTORY}/template_t1_denoised_gaussian.nii.gz
  --inputT2Volume ${ATLAS_INSTALL_DIRECTORY}/template_t2_denoised_gaussian.nii.gz
  --inputDiscreteVolume ${ATLAS_INSTALL_DIRECTORY}/template_ABC_labels.nii.gz
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/BRAINSContinuousClassTest_test.nii.gz
  )
