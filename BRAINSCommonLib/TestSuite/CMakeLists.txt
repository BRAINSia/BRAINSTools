add_executable(DumpImageInfo DumpImageInfo.cxx)
target_link_libraries(DumpImageInfo PRIVATE BRAINSCommonLib)
set_target_properties(DumpImageInfo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)
set_target_properties(DumpImageInfo PROPERTIES FOLDER ${MODULE_FOLDER})


# GTest-based unit tests for LocaleSafeConversions.h (issue #403)
# Uses BRAINSToolsGTest.cmake to locate/fetch GTest without SuperBuild changes.
include(BRAINSToolsGTest)
brainstools_enable_gtest()
add_executable(LocaleSafeConversionsGTest LocaleSafeConversionsGTest.cxx)
target_link_libraries(LocaleSafeConversionsGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(LocaleSafeConversionsGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
include(GoogleTest)
gtest_discover_tests(LocaleSafeConversionsGTest DISCOVERY_TIMEOUT 120)

# GTest port of PrettyPrintTableTest (issue #403)
add_executable(PrettyPrintTableGTest PrettyPrintTableGTest.cxx)
target_link_libraries(PrettyPrintTableGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(PrettyPrintTableGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(PrettyPrintTableGTest DISCOVERY_TIMEOUT 120)

# GTest port of Slicer3LandmarkIOExceptionTest (issue #403)
add_executable(Slicer3LandmarkIOGTest Slicer3LandmarkIOGTest.cxx)
target_link_libraries(Slicer3LandmarkIOGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(Slicer3LandmarkIOGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(Slicer3LandmarkIOGTest DISCOVERY_TIMEOUT 120)

# GTest port of AverageImageFilterTest (issue #403)
add_executable(AverageImageFilterGTest AverageImageFilterGTest.cxx)
target_link_libraries(AverageImageFilterGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(AverageImageFilterGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(AverageImageFilterGTest DISCOVERY_TIMEOUT 120)

# GTest coverage for DWIMetaDataDictionaryValidator (issue #403)
add_executable(DWIMetaDataDictionaryValidatorGTest DWIMetaDataDictionaryValidatorGTest.cxx)
target_link_libraries(DWIMetaDataDictionaryValidatorGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(DWIMetaDataDictionaryValidatorGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(DWIMetaDataDictionaryValidatorGTest DISCOVERY_TIMEOUT 120)

# GTest coverage for itkIO.h utilities (issue #403)
add_executable(itkIOGTest itkIOGTest.cxx)
target_link_libraries(itkIOGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(itkIOGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(itkIOGTest DISCOVERY_TIMEOUT 120)

# GTest port of ConvertToRigidAffine tests (issue #403)
add_executable(ConvertToRigidAffineGTest ConvertToRigidAffineGTest.cxx)
target_link_libraries(ConvertToRigidAffineGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(ConvertToRigidAffineGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(ConvertToRigidAffineGTest DISCOVERY_TIMEOUT 120)

# GTest coverage for Imgmath.h template functions (issue #403)
add_executable(ImgmathGTest ImgmathGTest.cxx)
target_link_libraries(ImgmathGTest PRIVATE BRAINSCommonLib GTest::gtest_main)
set_target_properties(ImgmathGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(ImgmathGTest DISCOVERY_TIMEOUT 120)

add_executable( itkResampleInPlaceImageFilterTest itkResampleInPlaceImageFilterTest.cxx)
set_target_properties(itkResampleInPlaceImageFilterTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)
target_link_libraries( itkResampleInPlaceImageFilterTest PRIVATE BRAINSCommonLib)
set_target_properties(itkResampleInPlaceImageFilterTest PROPERTIES FOLDER ${MODULE_FOLDER})

add_executable(BRAINSCleanMask BRAINSCleanMask.cxx)
target_link_libraries(BRAINSCleanMask PRIVATE BRAINSCommonLib)
set_target_properties(BRAINSCleanMask PROPERTIES FOLDER ${MODULE_FOLDER})

add_executable(DWIMetaDataDictionaryValidatorTests DWIMetaDataDictionaryValidatorTests.cxx)
target_link_libraries(DWIMetaDataDictionaryValidatorTests PRIVATE BRAINSCommonLib)
set_target_properties(DWIMetaDataDictionaryValidatorTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)
set_target_properties(DWIMetaDataDictionaryValidatorTests PROPERTIES FOLDER ${MODULE_FOLDER})

StandardBRAINSBuildMacro(NAME FindCenterOfBrain TARGET_LIBRARIES BRAINSCommonLib)

ExternalData_add_test(${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET}
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
ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME ${MetaDataValidatorTestName}
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

#if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 8) # This should be restored after fixing.
#ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET}
#  NAME FindCenterOfBrainOutputTest
#  COMMAND $<TARGET_FILE:ImageCompare>
#  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterDistanceImage.nii.gz
#  ${FindCenterDistanceImage}
#  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterGridImage.nii.gz
#  ${FindCenterGridImage}
#  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterAfterGridComputationsForegroundImage.nii.gz
#  ${FindCenterAfterGridComputationsForegroundImage}
#  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterClippedImageMask.nii.gz
#  ${FindCenterClippedImageMask}
#  ${CMAKE_CURRENT_BINARY_DIR}/FindCenterTrimmedImage.nii.gz
#  ${FindCenterTrimmedImage}
#  )
#set_property(TEST FindCenterOfBrainOutputTest APPEND PROPERTY DEPENDS
#  FindCenterOfBrainTest)
#
#
### HACK: INFO: Make this test do something meaningful.
### HACK: DAVID: This test is failing, but needs to have valid data added so that it will pass
#ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET}
#  NAME itkResampleInPlaceImageFilterTest
#  COMMAND $<TARGET_FILE:itkResampleInPlaceImageFilterTest>
#    itkResampleInPlaceImageFilterTest input1 transform1 checkresult
#)
#endif()
