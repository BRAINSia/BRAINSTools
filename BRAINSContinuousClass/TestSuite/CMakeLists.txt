###  This is the test suite for the BRAINSResample program.

MakeTestDriverFromSEMTool(BRAINSContinuousClass BRAINSContinuousClassTest.cxx)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSContinousClassTest_nii
  COMMAND ${LAUNCH_EXE} ${CMAKE_BIN_DIR}/BRAINSContinuousClassTestDriver
  --compare
  DATA{${TestData_DIR}/BRAINSContinuousClassTest_baseline.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSContinuousClassTest_test.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSContinuousClassTest
  --inputT1Volume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/template_t1.nii.gz
  --inputT2Volume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/template_t2.nii.gz
  --inputDiscreteVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/template_ABC_lables.nii.gz
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/BRAINSContinuousClassTest_test.nii.gz
  )

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
