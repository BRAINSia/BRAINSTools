configure_file( BRAINSABCSmallExtendedAtlasDefinition.xml.in
                ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallExtendedAtlasDefinition.xml
                @ONLY IMMEDIATE)


#add_executable(TestLinearRegressionTesting TestLinearRegresionInstensityMatching.cxx)
#target_link_libraries(TestLinearRegressionTesting BRAINSCommonLib )
#set_target_properties(TestLinearRegressionTesting PROPERTIES FOLDER ${MODULE_FOLDER})

MakeTestDriverFromSEMTool(BRAINSABC BRAINSABCTest.cxx)
target_include_directories(BRAINSABCTestDriver PRIVATE
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/brainseg
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/common
  ${BRAINSTools_BINARY_DIR}/BRAINSABC/brainseg
)
add_dependencies(BRAINSABCTestDriver InstallReferenceAtlas) ## Needed to ensure data is installed

# _brainstools_test_nproc is set once in BRAINSTools.cmake for all modules.

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME BRAINSABCSmallTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSABCTestDriver>
  --compare DATA{${TestData_DIR}/BRAINSABCSmallLabels.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallLabels.test.nii.gz
  --compareIntensityTolerance 1
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 10000
  BRAINSABCTest
   --atlasDefinition ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallExtendedAtlasDefinition.xml
   --atlasToSubjectInitialTransform DATA{${TestData_DIR}/BRAINSABCSmall_atlas_to_subject_transform.h5}
   --atlasToSubjectTransform BRAINSABCSmall_atlas_to_subject_transform.h5
   --atlasToSubjectTransformType Affine
   --debuglevel 0
   --filterIteration 0
   --filterMethod GradientAnisotropicDiffusion
   --gridSize 10,10,10
   --inputVolumeTypes T1,T2
   --inputVolumes DATA{${TestData_DIR}/affine_t1.nrrd}
   --inputVolumes DATA{${TestData_DIR}/affine_t2.nrrd}
   --interpolationMode Linear
   --maxBiasDegree 4
   --maxIterations 1
   --outputDir ./
   --outputDirtyLabels ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallvolume_label_seg.nii.gz
   --outputFormat NIFTI
   --outputLabels ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallLabels.test.nii.gz
   --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallT1_1.nii.gz
   --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallT2_1.nii.gz
   --posteriorTemplate ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallPOST_%s.nii.gz
   --purePlugsThreshold 0.2
)
set_tests_properties(BRAINSABCSmallTest PROPERTIES PROCESSORS ${_brainstools_test_nproc})

# Smoke variant on inputs shrunk by 2 (128^3); 64^3 is below BRAINSABC's minimum samples per label.
add_executable(BRAINSABCSmokeTestHelper BRAINSABCSmokeTestHelper.cxx)
target_link_libraries(BRAINSABCSmokeTestHelper PRIVATE ${BRAINSABC_ITK_LIBRARIES})
set_target_properties(BRAINSABCSmokeTestHelper PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})

foreach(_modality T1 T2)
  string(TOLOWER ${_modality} _modality_lower)
  ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME BRAINSABCSmokeShrink${_modality}Test
    COMMAND $<TARGET_FILE:BRAINSABCSmokeTestHelper> shrink 2
      DATA{${TestData_DIR}/affine_${_modality_lower}.nrrd}
      ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmoke${_modality}.nrrd
  )
  set_tests_properties(BRAINSABCSmokeShrink${_modality}Test PROPERTIES FIXTURES_SETUP BRAINSABCSmokeInputs)
endforeach()

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME BRAINSABCSmokeTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSABCTestDriver>
  BRAINSABCTest
   --atlasDefinition ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmallExtendedAtlasDefinition.xml
   --atlasToSubjectInitialTransform DATA{${TestData_DIR}/BRAINSABCSmall_atlas_to_subject_transform.h5}
   --atlasToSubjectTransform BRAINSABCSmoke_atlas_to_subject_transform.h5
   --atlasToSubjectTransformType Affine
   --debuglevel 0
   --filterIteration 0
   --filterMethod GradientAnisotropicDiffusion
   --gridSize 10,10,10
   --inputVolumeTypes T1,T2
   --inputVolumes ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokeT1.nrrd
   --inputVolumes ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokeT2.nrrd
   --interpolationMode Linear
   --maxBiasDegree 4
   --maxIterations 1
   --outputDir ./
   --outputDirtyLabels ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokeDirtyLabels.nii.gz
   --outputFormat NIFTI
   --outputLabels ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokeLabels.nii.gz
   --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokeT1_1.nii.gz
   --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokeT2_1.nii.gz
   --posteriorTemplate ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokePOST_%s.nii.gz
   --purePlugsThreshold 0.2
)
set_tests_properties(BRAINSABCSmokeTest PROPERTIES
  PROCESSORS ${_brainstools_test_nproc}
  FIXTURES_REQUIRED BRAINSABCSmokeInputs
  FIXTURES_SETUP BRAINSABCSmokeOutputs)

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME BRAINSABCSmokeLabelsTest
  COMMAND $<TARGET_FILE:BRAINSABCSmokeTestHelper> check
    DATA{${TestData_DIR}/BRAINSABCSmallLabels.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/BRAINSABCSmokeLabels.nii.gz 2 2
)
set_tests_properties(BRAINSABCSmokeLabelsTest PROPERTIES FIXTURES_REQUIRED BRAINSABCSmokeOutputs)

#if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 5) #These test takes way to long to run all the time
#ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME BRAINSABCLongTest
#  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSABCTestDriver>
#  --compare DATA{${TestData_DIR}/labels.nii.gz}
#  ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
#  --compareIntensityTolerance 1
#  --compareRadiusTolerance 1
#  --compareNumberOfPixelsTolerance 10000
#  BRAINSABCTest
#  --inputVolumes DATA{${TestData_DIR}/small_ISO_T1_REP0.nii.gz}
#  --inputVolumes DATA{${TestData_DIR}/small_ISO_T2_REP0.nii.gz}
#  --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/T1_cleaned.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T2_cleaned.nii.gz
#  --outputLabels ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
#  --outputDirtyLabels ${CMAKE_CURRENT_BINARY_DIR}/dirty_labels.nii.gz
#  --posteriorTemplate ${CMAKE_CURRENT_BINARY_DIR}/POST_%s.nii.gz
#  --inputVolumeTypes T1,T2
#  --filterIteration 3
#  --maxIterations 2
#  --maxBiasDegree 2
#  --debuglevel 0
#  --outputFormat NIFTI
#  --outputDir ${CMAKE_CURRENT_BINARY_DIR}
#  --gridSize  5,3,4
#  --atlasDefinition ${ATLAS_INSTALL_DIRECTORY}/ExtendedAtlasDefinition.xml
#  )
#ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME BRAINSABCLongTest2
#  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSABCTestDriver>
#  --compare DATA{${TestData_DIR}/labels.nii.gz}
#  ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
#  --compareIntensityTolerance 1
#  --compareRadiusTolerance 1
#  --compareNumberOfPixelsTolerance 10000
#  BRAINSABCTest
#  --inputVolumes DATA{${TestData_DIR}/small_ISO_T1_REP0.nii.gz}
#  --inputVolumes DATA{${TestData_DIR}/small_ISO_T1_REP0.nii.gz}
#  --inputVolumes DATA{${TestData_DIR}/small_ISO_T2_REP0.nii.gz}
#  --inputVolumes DATA{${TestData_DIR}/small_ISO_T2_REP0.nii.gz}
#  --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/T1_cleaned.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T1_cleaned2.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T2_cleaned.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T2_cleaned2.nii.gz
#  --outputLabels ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
#  --outputDirtyLabels ${CMAKE_CURRENT_BINARY_DIR}/dirty_labels.nii.gz
#  --posteriorTemplate ${CMAKE_CURRENT_BINARY_DIR}/POST_%s.nii.gz
#  --inputVolumeTypes T1,T1,T2,T2
#  --filterIteration 3
#  --maxIterations 2
#  --maxBiasDegree 2
#  --debuglevel 0
#  --outputFormat NIFTI
#  --outputDir ${CMAKE_CURRENT_BINARY_DIR}
#  --gridSize  5,3,4
#  --atlasDefinition ${ATLAS_INSTALL_DIRECTORY}/ExtendedAtlasDefinition.xml
#  )
#set_tests_properties(BRAINSABCLongTest PROPERTIES TIMEOUT 6500)
#endif()
#
# GTest port of BlendImageFilterTest (issue #403).
# The original was guarded behind BRAINSTools_MAX_TEST_LEVEL > 8; GTest makes
# the test unconditional and reports each failing pixel individually.
include(BRAINSToolsGTest)
brainstools_enable_gtest()
add_executable(BlendImageFilterGTest BlendImageFilterGTest.cxx)
target_include_directories(BlendImageFilterGTest PRIVATE
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/brainseg
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/common
  ${BRAINSTools_BINARY_DIR}/BRAINSABC/brainseg
)
target_link_libraries(BlendImageFilterGTest PRIVATE ${BRAINSABC_ITK_LIBRARIES} GTest::gtest_main)
set_target_properties(BlendImageFilterGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
include(GoogleTest)
gtest_discover_tests(BlendImageFilterGTest DISCOVERY_TIMEOUT 120)

add_executable(PosteriorFilenameTemplateGTest PosteriorFilenameTemplateGTest.cxx)
target_include_directories(PosteriorFilenameTemplateGTest PRIVATE ${BRAINSTools_SOURCE_DIR}/BRAINSABC/brainseg)
target_link_libraries(PosteriorFilenameTemplateGTest PRIVATE ${BRAINSABC_ITK_LIBRARIES} GTest::gtest_main)
set_target_properties(PosteriorFilenameTemplateGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(PosteriorFilenameTemplateGTest DISCOVERY_TIMEOUT 120)

add_executable(AtlasDefinitionGTest AtlasDefinitionGTest.cxx)
target_include_directories(AtlasDefinitionGTest PRIVATE ${BRAINSTools_SOURCE_DIR}/BRAINSABC/brainseg)
target_link_libraries(AtlasDefinitionGTest PRIVATE BRAINSABCCOMMONLIB GTest::gtest_main)
set_target_properties(AtlasDefinitionGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(AtlasDefinitionGTest DISCOVERY_TIMEOUT 120)

add_executable(OutputVolumeTemplateGTest OutputVolumeTemplateGTest.cxx)
target_include_directories(OutputVolumeTemplateGTest PRIVATE ${BRAINSTools_SOURCE_DIR}/BRAINSABC/brainseg)
target_link_libraries(OutputVolumeTemplateGTest PRIVATE ${BRAINSABC_ITK_LIBRARIES} GTest::gtest_main)
set_target_properties(OutputVolumeTemplateGTest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin
  FOLDER ${MODULE_FOLDER})
gtest_discover_tests(OutputVolumeTemplateGTest DISCOVERY_TIMEOUT 120)
