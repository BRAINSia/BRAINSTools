# --------------------------------------------------------- #
# BRAINSCut Testing
# --------------------------------------------------------- #
##
## BRAINSCut NetConfiguration(XML) file for testing
##
set( InputNetConfigurationForTesting ${BRAINSCut_SOURCE_DIR}/TestSuite/NetConfigurations/BRAINSCutTest.xml.in)
set( OutputNetConfigurationForTesting ${CMAKE_CURRENT_BINARY_DIR}/TestSuite/NetConfigurations/BRAINSCutTest.xml )

configure_file( ${InputNetConfigurationForTesting} ${OutputNetConfigurationForTesting}  IMMEDIATE @ONLY)

message("TestData_DIR=${TestData_DIR}")

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME FakeTestToDownloadFiles
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/FakeTest.cmake
  DATA{${TestData_DIR}/AtlasToSubjectScan1.${XFRM_EXT}}
  DATA{${TestData_DIR}/scan1T1.nii.gz}
  DATA{${TestData_DIR}/scan3T2.nii.gz}
  DATA{${TestData_DIR}/AtlasToSubjectScan2.${XFRM_EXT}}
  DATA{${TestData_DIR}/scan1T2.nii.gz}
  DATA{${TestData_DIR}/scan4Mask1.nii.gz}
  DATA{${TestData_DIR}/AtlasToSubjectScan3.${XFRM_EXT}}
  DATA{${TestData_DIR}/scan2Mask1.nii.gz}
  DATA{${TestData_DIR}/scan4T1.nii.gz}
  DATA{${TestData_DIR}/AtlasToSubjectScan4.${XFRM_EXT}}
  DATA{${TestData_DIR}/scan2T1.nii.gz}
  DATA{${TestData_DIR}/scan4T2.nii.gz}
  DATA{${TestData_DIR}/SubjectToAtlasScan1.${XFRM_EXT}}
  DATA{${TestData_DIR}/scan2T2.nii.gz}
  DATA{${TestData_DIR}/template_t1.nii.gz}
  DATA{${TestData_DIR}/SubjectToAtlasScan2.${XFRM_EXT}}
  DATA{${TestData_DIR}/scan3Mask1.nii.gz}
  DATA{${TestData_DIR}/scan1Mask1.nii.gz}
  DATA{${TestData_DIR}/scan3T1.nii.gz}
  )

##
## BRAINSCut Test Run
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ANNGenerateProbabilityTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --generateProbability --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ANNCreateVectorsTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --createVectors --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ANNTrainModelTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --trainModel --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ANNApplyModelTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --applyModel --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ANNAllInOneTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --generateProbability  --createVectors  --trainModel  --applyModel
  --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME RandomForestTrainModelTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --trainModel --method RandomForest 
  --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME RandomForestApplyModelTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --applyModel --method RandomForest
  --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME RandomForestTrainModelTestFromCLP
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --trainModel --method RandomForest --randomTreeDepth 6 --numberOfTrees 6
  --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME RandomForestApplyModelTestFromCLP
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --applyModel --method RandomForest --randomTreeDepth 6 --numberOfTrees 6
  --netConfiguration ${OutputNetConfigurationForTesting})

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

