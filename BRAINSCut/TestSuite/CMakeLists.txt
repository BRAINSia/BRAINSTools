# --------------------------------------------------------- #
# BRAINSCut Testing
# --------------------------------------------------------- #


##
## BRAINSCut NetConfiguration(XML) file for testing
##

add_executable(TestHashKey TestHashKey.cxx)
target_link_libraries(TestHashKey BRAINSCutCOMMONLIB)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME TestHashKeyUnitTests
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:TestHashKey> )

## ExternalData_expand_arguments( name variable_name_to_be_used file_downloaded?)

ExternalData_expand_arguments( ${PROJECT_NAME}FetchData AtlasToSubjectScan1 DATA{${TestData_DIR}/Transforms_h5/AtlasToSubjectScan1.${XFRM_EXT}} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan1T1     DATA{${TestData_DIR}/scan1T1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan3T2     DATA{${TestData_DIR}/scan3T2.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData AtlasToSubjectScan2 DATA{${TestData_DIR}/Transforms_h5/AtlasToSubjectScan2.${XFRM_EXT}} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan1T2     DATA{${TestData_DIR}/scan1T2.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan4Mask1  DATA{${TestData_DIR}/scan4Mask1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData AtlasToSubjectScan3 DATA{${TestData_DIR}/Transforms_h5/AtlasToSubjectScan3.${XFRM_EXT}} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan2Mask1  DATA{${TestData_DIR}/scan2Mask1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan4T1     DATA{${TestData_DIR}/scan4T1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData AtlasToSubjectScan4 DATA{${TestData_DIR}/Transforms_h5/AtlasToSubjectScan4.${XFRM_EXT}} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan2T1     DATA{${TestData_DIR}/scan2T1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan4T2     DATA{${TestData_DIR}/scan4T2.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData SubjectToAtlasScan1 DATA{${TestData_DIR}/Transforms_h5/SubjectToAtlasScan1.${XFRM_EXT}} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan2T2     DATA{${TestData_DIR}/scan2T2.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData template_t1 DATA{${TestData_DIR}/template_t1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData SubjectToAtlasScan2 DATA{${TestData_DIR}/Transforms_h5/SubjectToAtlasScan2.${XFRM_EXT}} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan3Mask1  DATA{${TestData_DIR}/scan3Mask1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan1Mask1  DATA{${TestData_DIR}/scan1Mask1.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan3T1     DATA{${TestData_DIR}/scan3T1.nii.gz} )

##
set( InputNetConfigurationForTesting ${BRAINSTools_SOURCE_DIR}/BRAINSCut/TestSuite/NetConfigurations/BRAINSCutTest.xml.in)
set( OutputNetConfigurationForTesting ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCutTest.xml )

configure_file( ${InputNetConfigurationForTesting} ${OutputNetConfigurationForTesting}  IMMEDIATE @ONLY)

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



##
## ordering of tests

## data dependency
set_property( TEST ANNGenerateProbabilityTest APPEND PROPERTY DEPENDS FakeTestToDownloadFiles)

## series of ANN test 
set_property( TEST ANNCreateVectorsTest       APPEND PROPERTY DEPENDS ANNGenerateProbabilityTest)
set_property( TEST ANNTrainModelTest          APPEND PROPERTY DEPENDS ANNCreateVectorsTest)
set_property( TEST ANNApplyModelTest          APPEND PROPERTY DEPENDS ANNTrainModelTest)

## data dependency
set_property( TEST ANNAllInOneTest            APPEND PROPERTY DEPENDS FakeTestToDownloadFiles)

## random forest
set_property( TEST RandomForestTrainModelTest  APPEND PROPERTY DEPENDS ANNCreateVectorsTest)
set_property( TEST RandomForestApplyModelTest  APPEND PROPERTY DEPENDS RandomForestTrainModelTest)
set_property( TEST RandomForestTrainModelTestFromCLP  APPEND PROPERTY DEPENDS ANNCreateVectorsTest)
set_property( TEST RandomForestApplyModelTestFromCLP  APPEND PROPERTY DEPENDS RandomForestTrainModelTestFromCLP)

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

add_subdirectory(src)
