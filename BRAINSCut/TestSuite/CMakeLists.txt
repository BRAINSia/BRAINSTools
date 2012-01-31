# --------------------------------------------------------- #
# BRAINSCut Testing
# --------------------------------------------------------- #
##
## BRAINSCut NetConfiguration(XML) file for testing
##
set( InputNetConfigurationForTesting ${BRAINSCut_SOURCE_DIR}/TestSuite/NetConfigurations/BRAINSCutTest.xml.in)
set( OutputNetConfigurationForTesting ${CMAKE_CURRENT_BINARY_DIR}/TestSuite/NetConfigurations/BRAINSCutTest.xml )

configure_file( ${InputNetConfigurationForTesting} ${OutputNetConfigurationForTesting}  IMMEDIATE @ONLY)

set(TestData_DIR ${CMAKE_CURRENT_SOURCE_DIR}/TestData/Input)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME FakeTestToDownloadFiles
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/FakeTest.cmake
  DATA{${TestData_DIR}/AtlasToSubjectScan1.mat.md5}
  DATA{${TestData_DIR}/scan1T1.nii.gz.md5}
  DATA{${TestData_DIR}/scan3T2.nii.gz.md5}
  DATA{${TestData_DIR}/AtlasToSubjectScan2.mat.md5}
  DATA{${TestData_DIR}/scan1T2.nii.gz.md5}
  DATA{${TestData_DIR}/scan4Mask1.nii.gz.md5}
  DATA{${TestData_DIR}/AtlasToSubjectScan3.mat.md5}
  DATA{${TestData_DIR}/scan2Mask1.nii.gz.md5}
  DATA{${TestData_DIR}/scan4T1.nii.gz.md5}
  DATA{${TestData_DIR}/AtlasToSubjectScan4.mat.md5}
  DATA{${TestData_DIR}/scan2T1.nii.gz.md5}
  DATA{${TestData_DIR}/scan4T2.nii.gz.md5}
  DATA{${TestData_DIR}/SubjectToAtlasScan1.mat.md5}
  DATA{${TestData_DIR}/scan2T2.nii.gz.md5}
  DATA{${TestData_DIR}/template_t1.nii.gz.md5}
  DATA{${TestData_DIR}/SubjectToAtlasScan2.mat.md5}
  DATA{${TestData_DIR}/scan3Mask1.nii.gz.md5}
  DATA{${TestData_DIR}/scan1Mask1.nii.gz.md5}
  DATA{${TestData_DIR}/scan3T1.nii.gz.md5}
  )

##
## BRAINSCut Test Run
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME GenerateProbabilityTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --generateProbability --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME CreateVectorsTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --createVectors --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME TrainModelTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --trainModel --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ApplyModelTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --applyModel --netConfiguration ${OutputNetConfigurationForTesting})
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME AllInOneTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --generateProbability  --createVectors  --trainModel  --applyModel
  --netConfiguration ${OutputNetConfigurationForTesting})

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target






