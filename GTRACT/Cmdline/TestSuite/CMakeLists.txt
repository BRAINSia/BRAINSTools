
#add_executable( gtractCoregBvaluesTests
#    ${GTRACT_SOURCE_DIR}/Cmdline/gtractCoregBvaluesPrimary.cxx gtractCoregBvaluesTests.cxx )
#target_link_libraries( gtractCoregBvaluesTests  BRAINSCommonLib GTRACTCommon )
#set_target_properties(gtractCoregBvaluesTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${GTRACT_BINARY_DIR})
#set(GTRACT_CoregBvalues_TESTS ${GTRACT_BINARY_DIR}/gtractCoregBvaluesTests )
#
#add_executable( gtractCoRegAnatomyTests
#    ${GTRACT_SOURCE_DIR}/Cmdline/gtractCoRegAnatomyPrimary.cxx gtractCoRegAnatomyTests.cxx )
#target_link_libraries( gtractCoRegAnatomyTests  BRAINSCommonLib GTRACTCommon )
#set_target_properties(gtractCoRegAnatomyTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${GTRACT_BINARY_DIR})
#set(GTRACT_CoRegAnatomy_TESTS ${GTRACT_BINARY_DIR}/gtractCoRegAnatomyTests )
#
#
#set(GTRACT_TEST_OUTPUT_DIR "${GTRACT_BINARY_DIR}/Testing")
#make_directory(${GTRACT_TEST_OUTPUT_DIR})
#
#set(GTRACT_TEST_INPUT_DIR "${GTRACT_SOURCE_DIR}/Data")
#
#message(warning "${CMAKE_BUNDLE_OUTPUT_DIRECTORY}")
#
#add_test(NAME gtractCoregBvalues_Test COMMAND ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoregBvalues_TESTS>
# --compare
#  ${GTRACT_SOURCE_DIR}/Data/b0_s.nhdr
#  ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
# --compareNumberOfPixelsTolerance 1000
# --compareIntensityTolerance 20
# gtractCoregBvaluesTest
#  --fixedVolume ${GTRACT_TEST_INPUT_DIR}/dti1.nhdr
#  --movingVolume ${GTRACT_TEST_INPUT_DIR}/dti2.nhdr
#  --fixedVolumeIndex 0 --debugLevel 0
#  --eddyCurrentCorrection --relaxationFactor 0.25 --maximumStepSize 0.1
#  --outputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
#  --outputTransform ${GTRACT_TEST_OUTPUT_DIR}/b0_s.${XFRM_EXT}
#)
#
#add_test(NAME gtractCoRegAnatomy_Rigid_Test COMMAND ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoRegAnatomy_TESTS>
#   gtractCoRegAnatomyTest
#     --inputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
#     --vectorIndex 0
#     --inputAnatomicalVolume ${GTRACT_TEST_INPUT_DIR}/t1_002.nhdr
#     --outputTransformName ${GTRACT_TEST_OUTPUT_DIR}/rigid_s.${XFRM_EXT}
#     --numberOfSamples 500000
#     --transformType Rigid
#)
#
#add_test(NAME gtractCoRegAnatomy_Bspline_Test COMMAND ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoRegAnatomy_TESTS>
#   gtractCoRegAnatomyTest
#     --inputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
#     --vectorIndex 0
#     --inputAnatomicalVolume ${GTRACT_TEST_INPUT_DIR}/t1_002.nhdr
#     --outputTransformName ${GTRACT_TEST_OUTPUT_DIR}/bspline_s.${XFRM_EXT}
#     --numberOfIterations 500
#     --transformType Bspline
#     --gradientTolerance 5e-6
#     --inputRigidTransform ${GTRACT_TEST_OUTPUT_DIR}/rigid_s.${XFRM_EXT}
#     --gridSize 18,6,18
#     --spatialScale 10
#)

add_executable( gtractResampleDWIInPlaceTests gtractResampleDWIInPlaceTests.cxx )
target_link_libraries( gtractResampleDWIInPlaceTests  BRAINSCommonLib GTRACTCommon )
set_target_properties(gtractResampleDWIInPlaceTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(GTRACTTestName GTRACTTest_gtractResampleDWIInPlace)
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${GTRACTTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:gtractResampleDWIInPlaceTests>
    --compare DATA{${TestData_DIR}/${GTRACTTestName}.standard.nrrd}
              ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
  gtractResampleDWIInPlaceTest
    --inputVolume DATA{${TestData_DIR}/DWI_test1.nrrd}
    --inputTransform DATA{${TestData_DIR}/${GTRACTTestName}_rigid_transform.mat}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${GTRACTTestName}.test.nrrd
)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )
