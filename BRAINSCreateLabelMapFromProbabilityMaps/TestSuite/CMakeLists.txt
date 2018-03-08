###  This is the test suite

MakeTestDriverFromSEMTool(BRAINSCreateLabelMapFromProbabilityMaps BRAINSCreateLabelMapFromProbabilityMapsTest.cxx)
add_dependencies(BRAINSCreateLabelMapFromProbabilityMapsTestDriver InstallReferenceAtlas) ## Needed to ensure data is installed

set(DIRTYLABELS ${CMAKE_CURRENT_BINARY_DIR}/DirtyLabels.nii.gz)
set(CLEANLABELS ${CMAKE_CURRENT_BINARY_DIR}/CleanLabels.nii.gz)

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET}
  NAME BRAINSCreateLabelMapFromProbabilityMapsTest
  COMMAND $<TARGET_FILE:BRAINSCreateLabelMapFromProbabilityMapsTestDriver>
  --compare ${DIRTYLABELS} DATA{${TestData_DIR}/BRAINSCreateLabelMapFromProbabilityMaps/DirtyLabels.nii.gz}
  --compare ${CLEANLABELS} DATA{${TestData_DIR}/BRAINSCreateLabelMapFromProbabilityMaps/CleanLabels.nii.gz}
  BRAINSCreateLabelMapFromProbabilityMapsTest
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_WM.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_SURFGM.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_BASAL.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_GLOBUS.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_THALAMUS.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_HIPPOCAMPUS.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_CRBLGM.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_CRBLWM.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_CSF.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_VB.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_NOTCSF.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_NOTGM.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_NOTWM.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_NOTVB.nii.gz
  --inputProbabilityVolume ${ATLAS_INSTALL_DIRECTORY}/GENERATED_AIR.nii.gz
  --priorLabelCodes  1,2,21,23,24,25,11,12,4,5,6,7,8,9,0
  --foregroundPriors 1,1,1,1,1,1,1,1,1,1,0,0,0,0,0
  --inclusionThreshold 0.0
  --dirtyLabelVolume ${DIRTYLABELS}
  --cleanLabelVolume ${CLEANLABELS}
)
