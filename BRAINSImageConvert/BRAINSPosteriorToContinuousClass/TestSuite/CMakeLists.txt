include_directories(
  ${BRAINSImageConvert_SOURCE_DIR}/BRAINSPosteriorToContinuousClass
  )

MakeTestDriverFromSEMTool(BRAINSPosteriorToContinuousClass BRAINSPosteriorToContinuousClassTest.cxx)

set(BRAINSImageConvertTestName BRAINSPosteriorToContinuousClassTest)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSImageConvertTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSImageConvertTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSImageConvertTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSImageConvertTestName}.output.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSPosteriorToContinuousClass
  --inputCrblWmVolume DATA{${ATLAS_SOURCE_DIR}/ALLPVCRBLWM.nii.gz}
  --inputCrblGmVolume DATA{${ATLAS_SOURCE_DIR}/ALLPVCRBLGM.nii.gz}
  --inputVbVolume DATA{${ATLAS_SOURCE_DIR}/ALLPVVB.nii.gz}
  --inputCsfVolume DATA{${ATLAS_SOURCE_DIR}/ALLPVCSF.nii.gz}
  --inputSurfaceGmVolume DATA{${ATLAS_SOURCE_DIR}/ALLPVSURFGM.nii.gz}
  --inputWhiteVolume DATA{${ATLAS_SOURCE_DIR}/ALLPVWM.nii.gz}
  --inputBasalGmVolume DATA{${ATLAS_SOURCE_DIR}/ALLPVBASALTISSUE.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSImageConvertTestName}.output.nii.gz}
#  --debugLevel 50
)