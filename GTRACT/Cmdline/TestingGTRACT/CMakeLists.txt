

add_executable( gtractCoregBvaluesTests
    ${GTRACT_SOURCE_DIR}/Cmdline/gtractCoregBvaluesPrimary.cxx gtractCoregBvaluesTests.cxx )
target_link_libraries( gtractCoregBvaluesTests  BRAINSCommonLib GTRACTCommon )
set_target_properties(gtractCoregBvaluesTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${GTRACT_BINARY_DIR})
set(GTRACT_CoregBvalues_TESTS ${GTRACT_BINARY_DIR}/gtractCoregBvaluesTests )

add_executable( gtractCoRegAnatomyTests
    ${GTRACT_SOURCE_DIR}/Cmdline/gtractCoRegAnatomyPrimary.cxx gtractCoRegAnatomyTests.cxx )
target_link_libraries( gtractCoRegAnatomyTests  BRAINSCommonLib GTRACTCommon )
set_target_properties(gtractCoRegAnatomyTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${GTRACT_BINARY_DIR})
set(GTRACT_CoRegAnatomy_TESTS ${GTRACT_BINARY_DIR}/gtractCoRegAnatomyTests )


set(GTRACT_TEST_OUTPUT_DIR "${GTRACT_BINARY_DIR}/Testing")
make_directory(${GTRACT_TEST_OUTPUT_DIR})

set(GTRACT_TEST_INPUT_DIR "${GTRACT_SOURCE_DIR}/Data")

message(warning "${CMAKE_BUNDLE_OUTPUT_DIRECTORY}")

#add_test(gtractCoregBvalues_Test ${CMAKE_BUNDLE_OUTPUT_DIRECTORY}/gtractCoregBvalues
add_test(gtractCoregBvalues_Test ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoregBvalues_TESTS>
 --compare
  ${GTRACT_SOURCE_DIR}/Data/b0_s.nhdr
  ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
 --compareNumberOfPixelsTolerance 1000
 --compareIntensityTolerance 20
 gtractCoregBvaluesTest
  --fixedVolume ${GTRACT_TEST_INPUT_DIR}/dti1.nhdr
  --movingVolume ${GTRACT_TEST_INPUT_DIR}/dti2.nhdr
  --fixedVolumeIndex 0 --debugLevel 0
  --eddyCurrentCorrection --relaxationFactor 0.25 --maximumStepSize 0.1
  --outputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
  --outputTransform ${GTRACT_TEST_OUTPUT_DIR}/b0_s.mat
)

add_test(gtractCoRegAnatomy_Rigid_Test ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoRegAnatomy_TESTS>
   gtractCoRegAnatomyTest
     --inputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
     --vectorIndex 0
     --inputAnatomicalVolume ${GTRACT_TEST_INPUT_DIR}/t1_002.nhdr
     --outputTransformName ${GTRACT_TEST_OUTPUT_DIR}/rigid_s.mat
     --numberOfSamples 500000
     --transformType Rigid
)

add_test(gtractCoRegAnatomy_Bspline_Test ${LAUNCH_EXE} $<TARGET_FILE:GTRACT_CoRegAnatomy_TESTS>
   gtractCoRegAnatomyTest
     --inputVolume ${GTRACT_TEST_OUTPUT_DIR}/b0_s.nhdr
     --vectorIndex 0
     --inputAnatomicalVolume ${GTRACT_TEST_INPUT_DIR}/t1_002.nhdr
     --outputTransformName ${GTRACT_TEST_OUTPUT_DIR}/bspline_s.mat
     --numberOfIterations 500
     --transformType Bspline
     --gradientTolerance 5e-6
     --inputRigidTransform ${GTRACT_TEST_OUTPUT_DIR}/rigid_s.mat
     --gridSize 18,6,18
     --spatialScale 10
)

