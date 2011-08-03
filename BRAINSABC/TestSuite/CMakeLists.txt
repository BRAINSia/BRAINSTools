list(APPEND ExternalData_URL_TEMPLATES
  # Local data store populated by the ITK pre-commit hook
  "file:///${CMAKE_SOURCE_DIR}/.ExternalData/%(algo)/%(hash)"
  # Data published by Iowa Psychiatry web interface
  "http://www.psychiatry.uiowa.edu/users/brainstestdata/ctestdata/%(algo)/%(hash)"

  # Data published by MIDAS
  "http://midas.kitware.com/api/rest/midas.bitstream.by.hash?hash=%(hash)&algorithm=%(algo)"

  # Data published by developers using git-gerrit-push.
  "http://www.itk.org/files/ExternalData/%(algo)/%(hash)"
  )

# Tell ExternalData commands to transform raw files to content links.
# TODO: Condition this feature on presence of our pre-commit hook.
set(ExternalData_LINK_CONTENT MD5)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/ExternalData.cmake)

include_directories(
  ${BRAINSABC_SOURCE_DIR}/brainseg
  ${BRAINSABC_SOURCE_DIR}/common
  ${BRAINSABC_BINARY_DIR}/brainseg
)



MakeTestDriverFromSEMTool(BRAINSABC BRAINSABCTest.cxx)

ExternalData_add_test( FetchData NAME BRAINSABCLongTest
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSABCTestDriver>
    --compare DATA{test_data/refResults/labels.nii.gz}
               ${CMAKE_CURRENT_BINARY_DIR}/labels.nii.gz
    --compareIntensityTolerance 1
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 10000
  BRAINSABCTest
  --inputVolumes DATA{test_data/small_ISO_T1_REP0.nii.gz} --inputVolumes DATA{test_data/small_ISO_T2_REP0.nii.gz}
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
  --atlasDefinition ${CMAKE_CURRENT_BINARY_DIR}/../../src/bin/Atlas/Atlas_20110607/AtlasPVDefinition.xml
)
set_tests_properties(BRAINSABCLongTest PROPERTIES TIMEOUT 6500)

if(0) # This should be restored after fixing.
add_executable(BlendImageFilterTest BlendImageFilterTest.cxx)
target_link_libraries(BlendImageFilterTest ${ITK_LIBRARIES})
ExternalData_add_test( FetchData NAME BlendImageFilterTest COMMAND ${LAUNCH_EXE}  $<TARGET_FILE:BlendImageFilterTest> )
endif()

ExternalData_Add_Target( FetchData )  # Name of data management target
