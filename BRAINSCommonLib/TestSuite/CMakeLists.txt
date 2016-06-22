
add_executable(PrettyPrintTableTest PrettyPrintTableTest.cxx)
target_link_libraries(PrettyPrintTableTest BRAINSCommonLib)
set_target_properties(PrettyPrintTableTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)

add_executable( itkResampleInPlaceImageFilterTest itkResampleInPlaceImageFilterTest.cxx)
set_target_properties(itkResampleInPlaceImageFilterTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)
target_link_libraries( itkResampleInPlaceImageFilterTest ${BRAINSCommonLib_ITK_LIBRARIES})

add_executable(BRAINSCleanMask BRAINSCleanMask.cxx)
target_link_libraries(BRAINSCleanMask ${BRAINSCommonLib_ITK_LIBRARIES})

add_executable(DWIMetaDataDictionaryValidatorTests DWIMetaDataDictionaryValidatorTests.cxx)
target_link_libraries(DWIMetaDataDictionaryValidatorTests BRAINSCommonLib)
set_target_properties(DWIMetaDataDictionaryValidatorTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)

add_executable(Slicer3LandmarkIOExceptionTest Slicer3LandmarkIOTest.cxx)
target_link_libraries(Slicer3LandmarkIOExceptionTest BRAINSCommonLib)
set_target_properties(Slicer3LandmarkIOExceptionTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)
StandardBRAINSBuildMacro(NAME FindCenterOfBrain TARGET_LIBRARIES BRAINSCommonLib)

ExternalData_add_test(FindCenterOfBrainFetchData
  NAME PrettyPrintTableTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:PrettyPrintTableTest>
  ## No arguments
  )

ExternalData_add_test(FindCenterOfBrainFetchData
  NAME Slicer3LandmarkIOExceptionTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:Slicer3LandmarkIOExceptionTest>
  ## No arguments
  )

ExternalData_add_test(FindCenterOfBrainFetchData
  NAME FindCenterOfBrainTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:FindCenterOfBrain>
  --inputVolume DATA{${TestData_DIR}/BRAINSEyeDetectorTest_T1_standard.nii.gz}
  --generateDebugImages
  --debugDistanceImage ${CMAKE_CURRENT_BINARY_DIR}/FindCenterDistanceImage.nii.gz
  --debugGridImage ${CMAKE_CURRENT_BINARY_DIR}/FindCenterGridImage.nii.gz
  --debugAfterGridComputationsForegroundImage
  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterAfterGridComputationsForegroundImage.nii.gz
  --debugClippedImageMask ${CMAKE_CURRENT_BINARY_DIR}/FindCenterClippedImageMask.nii.gz
  --debugTrimmedImage ${CMAKE_CURRENT_BINARY_DIR}/FindCenterTrimmedImage.nii.gz
  )

ExternalData_expand_arguments(FindCenterOfBrainFetchData FindCenterDistanceImage
  DATA{${TestData_DIR}/FindCenter/FindCenterDistanceImage.nii.gz})
ExternalData_expand_arguments(FindCenterOfBrainFetchData FindCenterGridImage
  DATA{${TestData_DIR}/FindCenter/FindCenterGridImage.nii.gz})

ExternalData_expand_arguments(FindCenterOfBrainFetchData FindCenterAfterGridComputationsForegroundImage
  DATA{${TestData_DIR}/FindCenter/FindCenterAfterGridComputationsForegroundImage.nii.gz})

ExternalData_expand_arguments(FindCenterOFBrainFetchData FindCenterClippedImageMask
  DATA{${TestData_DIR}/FindCenter/FindCenterClippedImageMask.nii.gz})

ExternalData_expand_arguments(FindCenterOfBrainFetchData
  FindCenterTrimmedImage DATA{${TestData_DIR}/FindCenter/FindCenterTrimmedImage.nii.gz})

set(MetaDataValidatorTestName DWIMetaDataDictionaryValidatorTest)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${MetaDataValidatorTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:DWIMetaDataDictionaryValidatorTests>
    --compare DATA{${TestData_DIR}/DWI_TestData_OUTPUTS/PhilipsAchieva2.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/outputSameAsRefVolume.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  DWIMetaDataDictionaryValidatorTest
  ${CMAKE_CURRENT_BINARY_DIR}/ouputArtificialTestNrrdImage.nrrd
  DATA{${TestData_DIR}/DWI_TestData_OUTPUTS/PhilipsAchieva2.nrrd}
  ${CMAKE_CURRENT_BINARY_DIR}/outputSameAsRefVolume.nrrd
)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )

if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 8) # This should be restored after fixing.
ExternalData_add_test(FindCenterOfBrainFetchData
  NAME FindCenterOfBrainOutputTest
  COMMAND $<TARGET_FILE:ImageCompare>
  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterDistanceImage.nii.gz
  ${FindCenterDistanceImage}
  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterGridImage.nii.gz
  ${FindCenterGridImage}
  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterAfterGridComputationsForegroundImage.nii.gz
  ${FindCenterAfterGridComputationsForegroundImage}
  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterClippedImageMask.nii.gz
  ${FindCenterClippedImageMask}
  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterTrimmedImage.nii.gz
  ${FindCenterTrimmedImage}
  )
set_property(TEST FindCenterOfBrainOutputTest APPEND PROPERTY DEPENDS
  FindCenterOfBrainTest)


## HACK: TODO: Make this test do something meaningful.
## HACK: DAVID: This test is failing, but needs to have valid data added so that it will pass
ExternalData_add_test(FindCenterOfBrainFetchData
  NAME itkResampleInPlaceImageFilterTest
  COMMAND $<TARGET_FILE:itkResampleInPlaceImageFilterTest>
    itkResampleInPlaceImageFilterTest input1 transform1 checkresult
)
endif()
ExternalData_add_target(FindCenterOfBrainFetchData)
