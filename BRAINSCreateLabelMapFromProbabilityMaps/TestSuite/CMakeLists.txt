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
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_WM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_SURFGM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_ACCUMBEN.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_CAUDATE.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_PUTAMEN.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_GLOBUS.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_THALAMUS.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_HIPPOCAMPUS.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_CRBLGM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_CRBLWM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_CSF.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_VB.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_NOTCSF.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_NOTGM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_NOTWM.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_NOTVB.nii.gz
  --inputProbabilityVolume ${REF_DIR}/EXTENDED_AIR.nii.gz
  --priorLabelCodes  1,2,20,21,22,23,24,25,11,12,4,5,6,7,8,9,0
  --foregroundPriors 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0
  --inclusionThreshold 0.0
  --dirtyLabelVolume ${DIRTYLABELS}
  --cleanLabelVolume ${CLEANLABELS}
)
ExternalData_Add_Target( BRAINSCreateLabelMapFromProbabilityMapsFetchData )