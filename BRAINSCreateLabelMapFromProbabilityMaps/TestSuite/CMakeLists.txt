MakeTestDriverFromSEMTool(BRAINSCreateLabelMapFromProbabilityMaps
  BRAINSCreateLabelMapFromProbabilityMapsTest.cxx)

set(REF_DIR ${ReferenceAtlas_DIR}/${ATLAS_NAME})

set(DIRTYLABELS ${CMAKE_CURRENT_BINARY_DIR}/DirtyLabels.nii.gz)
set(CLEANLABELS ${CMAKE_CURRENT_BINARY_DIR}/CleanLabels.nii.gz)

ExternalData_add_test(BRAINSCreateLabelMapFromProbabilityMapsFetchData
  NAME BRAINSCreateLabelMapFromProbabilityMapsTest
  COMMAND $<TARGET_FILE:BRAINSCreateLabelMapFromProbabilityMapsTestDriver>
  --compare ${DIRTYLABELS} DATA{${TestData_DIR}/BRAINSCreateLabelMapFromProbabilityMaps/DirtyLabels.nii.gz}
  --compare ${CLEANLABELS} DATA{${TestData_DIR}/BRAINSCreateLabelMapFromProbabilityMaps/CleanLabels.nii.gz}
  BRAINSCreateLabelMapFromProbabilityMapsTest
  --inputProbabilityVolume ${REF_DIR}/GENERATED_WM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_SURFGM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_BASAL.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_GLOBUS.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_THALAMUS.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_HIPPOCAMPUS.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_CRBLGM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_CRBLWM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_CSF.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_VB.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_NOTCSF.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_NOTGM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_NOTWM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_NOTVB.nii.gz
  --inputProbabilityVolume ${REF_DIR}/GENERATED_AIR.nii.gz
  --priorLabelCodes  1,2,21,23,24,25,11,12,4,5,6,7,8,9,0
  --foregroundPriors 1,1,1,1,1,1,1,1,1,1,0,0,0,0,0
  --inclusionThreshold 0.0
  --dirtyLabelVolume ${DIRTYLABELS}
  --cleanLabelVolume ${CLEANLABELS}
)
ExternalData_Add_Target( BRAINSCreateLabelMapFromProbabilityMapsFetchData )
