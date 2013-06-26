include_directories(
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/brainseg
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/common
  ${BRAINSTools_BINARY_DIR}/BRAINSABC/brainseg
  )



MakeTestDriverFromSEMTool(BRAINSABC BRAINSABCTest.cxx)

if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 5) #These test takes way to long to run all the time
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSABCLongTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSABCTestDriver>
  --compare DATA{${TestData_DIR}/labels.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
  --compareIntensityTolerance 1
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 10000
  BRAINSABCTest
  --inputVolumes DATA{${TestData_DIR}/small_ISO_T1_REP0.nii.gz}
  --inputVolumes DATA{${TestData_DIR}/small_ISO_T2_REP0.nii.gz}
  --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/T1_cleaned.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T2_cleaned.nii.gz
  --outputLabels ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
  --outputDirtyLabels ${CMAKE_CURRENT_BINARY_DIR}/dirty_labels.nii.gz
  --posteriorTemplate ${CMAKE_CURRENT_BINARY_DIR}/POST_%s.nii.gz
  --inputVolumeTypes T1,T2
  --filterIteration 3
  --maxIterations 2
  --maxBiasDegree 2
  --debuglevel 0
  --outputFormat NIFTI
  --outputDir ${CMAKE_CURRENT_BINARY_DIR}
  --gridSize  5,3,4
  --atlasDefinition ${ReferenceAtlas_DIR}/${ATLAS_NAME}/AtlasPVDefinition.xml
  )
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSABCLongTest2
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSABCTestDriver>
  --compare DATA{${TestData_DIR}/labels.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
  --compareIntensityTolerance 1
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 10000
  BRAINSABCTest
  --inputVolumes DATA{${TestData_DIR}/small_ISO_T1_REP0.nii.gz}
  --inputVolumes DATA{${TestData_DIR}/small_ISO_T1_REP0.nii.gz}
  --inputVolumes DATA{${TestData_DIR}/small_ISO_T2_REP0.nii.gz}
  --inputVolumes DATA{${TestData_DIR}/small_ISO_T2_REP0.nii.gz}
  --outputVolumes ${CMAKE_CURRENT_BINARY_DIR}/T1_cleaned.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T1_cleaned2.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T2_cleaned.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/T2_cleaned2.nii.gz
  --outputLabels ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
  --outputDirtyLabels ${CMAKE_CURRENT_BINARY_DIR}/dirty_labels.nii.gz
  --posteriorTemplate ${CMAKE_CURRENT_BINARY_DIR}/POST_%s.nii.gz
  --inputVolumeTypes T1,T1,T2,T2
  --filterIteration 3
  --maxIterations 2
  --maxBiasDegree 2
  --debuglevel 0
  --outputFormat NIFTI
  --outputDir ${CMAKE_CURRENT_BINARY_DIR}
  --gridSize  5,3,4
  --atlasDefinition ${ReferenceAtlas_DIR}/${ATLAS_NAME}/AtlasPVDefinition.xml
  )
set_tests_properties(BRAINSABCLongTest PROPERTIES TIMEOUT 6500)
endif()

if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 8) # This should be restored after fixing.
  add_executable(BlendImageFilterTest BlendImageFilterTest.cxx)
  target_link_libraries(BlendImageFilterTest ${ITK_LIBRARIES})
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BlendImageFilterTest COMMAND ${LAUNCH_EXE}  $<TARGET_FILE:BlendImageFilterTest> )
endif()

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
