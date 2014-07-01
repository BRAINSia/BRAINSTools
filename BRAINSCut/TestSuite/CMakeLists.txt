# --------------------------------------------------------- #
# BRAINSCut Testing
# --------------------------------------------------------- #

##//Testing directory uses the ExternalData directory
set(Reference_BINARY_DIR ${BRAINSTools_BINARY_DIR}/ExternalData)

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
ExternalData_expand_arguments( ${PROJECT_NAME}FetchData scan3SG     DATA{${TestData_DIR}/scan3SG.nii.gz} )
ExternalData_expand_arguments( ${PROJECT_NAME}FetchData BRAINSCutTemplateToScan3     DATA{${TestData_DIR}/BRAINSCutTemplateToScan3.h5} )

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

## machine-learning model stability test
# - for a identical target image, a given deformation, and a given model,
#   the segmentation results should be identical.

# caudate
set( InputNetConfigurationForOpenCVTesting_caudate
  ${BRAINSTools_SOURCE_DIR}/BRAINSCut/TestSuite/NetConfigurations/OpenCVStabilityTest_caudate.xml.in)
set( OutputNetConfigurationForOpenCVTesting_caudate
  ${CMAKE_CURRENT_BINARY_DIR}/OpenCVStabilityTest_caudate.xml )

configure_file( ${InputNetConfigurationForOpenCVTesting_caudate} ${OutputNetConfigurationForOpenCVTesting_caudate}  IMMEDIATE @ONLY)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME OpenCVConsistencyTest_caudate
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --applyModel --method RandomForest
  --modelFilename ${Reference_BINARY_DIR}/modelFiles/trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz
  --netConfiguration ${OutputNetConfigurationForOpenCVTesting_caudate})
# putamen
set( InputNetConfigurationForOpenCVTesting_putamen
  ${BRAINSTools_SOURCE_DIR}/BRAINSCut/TestSuite/NetConfigurations/OpenCVStabilityTest_putamen.xml.in)
set( OutputNetConfigurationForOpenCVTesting_putamen
  ${CMAKE_CURRENT_BINARY_DIR}/OpenCVStabilityTest_putamen.xml )

configure_file( ${InputNetConfigurationForOpenCVTesting_putamen} ${OutputNetConfigurationForOpenCVTesting_putamen}  IMMEDIATE @ONLY)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME OpenCVConsistencyTest_putamen
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSCut>
  --applyModel --method RandomForest
  --modelFilename ${Reference_BINARY_DIR}/modelFiles/trainModelFile.txtD0060NT0060_putamen.gz
  --netConfiguration ${OutputNetConfigurationForOpenCVTesting_putamen})

# - As a number of trees increases, the training error should be decreased
#   approximately same rate to the previously trained one.
# TODO


##
## ordering of tests

## data dependency
set_property( TEST ANNGenerateProbabilityTest APPEND PROPERTY DEPENDS ${PROJECT_NAME}FetchData)

## series of ANN test
set_property( TEST ANNCreateVectorsTest       APPEND PROPERTY DEPENDS ANNGenerateProbabilityTest)
set_property( TEST ANNTrainModelTest          APPEND PROPERTY DEPENDS ANNCreateVectorsTest)
set_property( TEST ANNApplyModelTest          APPEND PROPERTY DEPENDS ANNTrainModelTest)

## data dependency
set_property( TEST ANNAllInOneTest            APPEND PROPERTY DEPENDS ANNApplyModelTest)

## random forest
set_property( TEST RandomForestTrainModelTest  APPEND PROPERTY DEPENDS ANNAllInOneTest)
set_property( TEST RandomForestApplyModelTest  APPEND PROPERTY DEPENDS RandomForestTrainModelTest)
set_property( TEST RandomForestTrainModelTestFromCLP  APPEND PROPERTY DEPENDS RandomForestApplyModelTest)
set_property( TEST RandomForestApplyModelTestFromCLP  APPEND PROPERTY DEPENDS RandomForestTrainModelTestFromCLP)

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

add_subdirectory(src)
