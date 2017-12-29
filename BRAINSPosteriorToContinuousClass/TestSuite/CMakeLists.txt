include_directories(
  ${BRAINSPosteriorToContinuousClass_SOURCE_DIR}/BRAINSPosteriorToContinuousClass
  )

MakeTestDriverFromSEMTool(BRAINSPosteriorToContinuousClass BRAINSPosteriorToContinuousClassTest.cxx)
add_dependencies(BRAINSPosteriorToContinuousClassTestDriver InstallReferenceAtlas) ## Needed to ensure data is installed

set(BRAINSPosteriorToContinuousClassTestName BRAINSPosteriorToContinuousClassTest)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSPosteriorToContinuousClassTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSPosteriorToContinuousClassTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSPosteriorToContinuousClassTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSPosteriorToContinuousClassTestName}.output.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSPosteriorToContinuousClassTest
  --inputCrblWmVolume ${ATLAS_INSTALL_DIRECTORY}/ALLPVCRBLWM.nii.gz
  --inputCrblGmVolume ${ATLAS_INSTALL_DIRECTORY}/ALLPVCRBLGM.nii.gz
  --inputVbVolume ${ATLAS_INSTALL_DIRECTORY}/ALLPVVB.nii.gz
  --inputCsfVolume ${ATLAS_INSTALL_DIRECTORY}/ALLPVCSF.nii.gz
  --inputSurfaceGmVolume ${ATLAS_INSTALL_DIRECTORY}/ALLPVSURFGM.nii.gz
  --inputWhiteVolume ${ATLAS_INSTALL_DIRECTORY}/ALLPVWM.nii.gz
  --inputBasalGmVolume ${ATLAS_INSTALL_DIRECTORY}/ALLPVBASALTISSUE.nii.gz
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSPosteriorToContinuousClassTestName}.output.nii.gz
#  --debugLevel 50
)
