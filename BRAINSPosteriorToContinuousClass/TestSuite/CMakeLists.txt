include_directories(
  ${BRAINSPosteriorToContinuousClass_SOURCE_DIR}/BRAINSPosteriorToContinuousClass
  )

MakeTestDriverFromSEMTool(BRAINSPosteriorToContinuousClass BRAINSPosteriorToContinuousClassTest.cxx)

set(BRAINSPosteriorToContinuousClassTestName BRAINSPosteriorToContinuousClassTest)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSPosteriorToContinuousClassTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSPosteriorToContinuousClass>
  --compare DATA{${TestData_DIR}/${BRAINSPosteriorToContinuousClassTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSPosteriorToContinuousClassTestName}.output.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSPosteriorToContinuousClass
  --inputCrblWmVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/ALLPVCRBLWM.nii.gz
  --inputCrblGmVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/ALLPVCRBLGM.nii.gz
  --inputVbVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/ALLPVVB.nii.gz
  --inputCsfVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/ALLPVCSF.nii.gz
  --inputSurfaceGmVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/ALLPVSURFGM.nii.gz
  --inputWhiteVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/ALLPVWM.nii.gz
  --inputBasalGmVolume ${ReferenceAtlas_DIR}/${ATLAS_NAME}/ALLPVBASALTISSUE.nii.gz
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSPosteriorToContinuousClassTestName}.output.nii.gz
#  --debugLevel 50
)
