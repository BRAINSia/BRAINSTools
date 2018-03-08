#
# ITK
#

FindITKUtil( BRAINSCut_TestSuite_ITK
  ITKImageCompare
  ITKTestKernel
)
## TESTING EX.
## [BRAINSTools_Build_Dir]/bin/testRandomForest [correct_dir]/optdigits.train [correct_dir]/optdigits.test [output_dir]/optdigitsTrainedModel.txt
# Below Commented out by Regina
# Below Commented out by Regina
# Below Commented out by Regina
# Below Commented out by Regina
#add_executable(testRandomForest testRandomForest.cxx)
#target_link_libraries(testRandomForest ${BRAINSCut_TestSuite_ITK_LIBRARIES} ${OpenCV_LIBRARIES})
#set_target_properties(testRandomForest PROPERTIES FOLDER ${MODULE_FOLDER})
#
#ExternalData_add_test( ${PROJECT_NAME}FetchData NAME OpenCVRandomForestTest
#  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:testRandomForest>
#  DATA{${TestData_DIR}/optdigits.train}
#  DATA{${TestData_DIR}/optdigits.test}
#  ${CMAKE_CURRENT_BINARY_DIR}/optdigitTrainedModel.output )

add_executable(compareBinaryTest compareBinaryTest.cxx)
target_link_libraries(compareBinaryTest ${BRAINSCut_TestSuite_ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES})
set_target_properties(compareBinaryTest PROPERTIES FOLDER ${MODULE_FOLDER})

# compare caudate
ExternalData_expand_arguments( ${PROJECT_NAME}FetchData BRAINSCutcaudate_left_scan3_baseline
    DATA{${TestData_DIR}/BRAINSCutcaudate_left_scan3_baseline.nii.gz} )
ExternalData_expand_arguments( ${PROJECT_NAME}FetchData BRAINSCutcaudate_right_scan3_baseline
    DATA{${TestData_DIR}/BRAINSCutcaudate_right_scan3_baseline.nii.gz} )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME VolumeMatch_Left_Caudate
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:compareBinaryTest>
    ${CMAKE_CURRENT_BINARY_DIR}/../l_caudate_scan3.nii.gz
    ${BRAINSCutcaudate_left_scan3_baseline} )
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME VolumeMatch_Right_Caudate
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:compareBinaryTest>
    ${CMAKE_CURRENT_BINARY_DIR}/../r_caudate_scan3.nii.gz
    ${BRAINSCutcaudate_right_scan3_baseline} )

# comare putamen
ExternalData_expand_arguments( ${PROJECT_NAME}FetchData BRAINSCutputamen_left_scan3_baseline
    DATA{${TestData_DIR}/BRAINSCutputamen_left_scan3_baseline.nii.gz} )
  ExternalData_expand_arguments( ${PROJECT_NAME}FetchData BRAINSCutputamen_right_scan3_baseline
    DATA{${TestData_DIR}/BRAINSCutputamen_right_scan3_baseline.nii.gz} )
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME VolumeMatch_Left_Putamen
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:compareBinaryTest>
    ${CMAKE_CURRENT_BINARY_DIR}/../l_putamen_scan3.nii.gz
    ${BRAINSCutputamen_left_scan3_baseline} )
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME VolumeMatch_Right_Putamen
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:compareBinaryTest>
    ${CMAKE_CURRENT_BINARY_DIR}/../r_putamen_scan3.nii.gz
    ${BRAINSCutputamen_right_scan3_baseline} )

  set_property(TEST VolumeMatch_Left_Caudate  APPEND PROPERTY DEPENDS BRAINSCutcaudate_left_scan3_baseline)
  set_property(TEST VolumeMatch_Left_Caudate  APPEND PROPERTY DEPENDS OpenCVConsistencyTest_caudate)
  set_property(TEST VolumeMatch_Right_Caudate APPEND PROPERTY DEPENDS BRAINSCutcaudate_right_scan3_baseline)
  set_property(TEST VolumeMatch_Right_Caudate APPEND PROPERTY DEPENDS OpenCVConsistencyTest_caudate)

  set_property(TEST VolumeMatch_Left_Putamen  APPEND PROPERTY DEPENDS BRAINSCutputamen_left_scan3_baseline)
  set_property(TEST VolumeMatch_Left_Putamen  APPEND PROPERTY DEPENDS OpenCVConsistencyTest_putamen)
  set_property(TEST VolumeMatch_Right_Putamen APPEND PROPERTY DEPENDS BRAINSCutputamen_right_scan3_baseline)
  set_property(TEST VolumeMatch_Right_Putamen APPEND PROPERTY DEPENDS OpenCVConsistencyTest_putamen)
