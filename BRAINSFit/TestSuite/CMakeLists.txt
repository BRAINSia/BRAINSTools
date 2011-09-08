
MakeTestDriverFromSEMTool(BRAINSFit BRAINSFitTest.cxx)

if(COMPILE_TEST_GENERATION_PROGRAMS)  ## These were programs needed to create the test suite, but are not needed for actual testing!
  # The following program was used to generate TestData/test2.nii.gz
  # from test.nii.gz
  if(MakeMakeXfrmImage)
    add_executable(makexfrmedImage makexfrmedImage.cxx)
    target_link_libraries(makexfrmedImage ${ITK_LIBRARIES})
  endif(MakeMakeXfrmImage)

  if(BRAINSFitTestDriver)
    # The test driver is only needed to help generate new testcases
    add_executable(BRAINSFitTestDriver BRAINSFitTestDriver.cxx)
    target_link_libraries(BRAINSFitTestDriver ${ITK_LIBRARIES}
      ${OPTIONAL_DEBUG_LINK_LIBRARIES} )
    set_target_properties(BRAINSFitTestDriver PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BRAINSFit_BINARY_DIR})
  endif()
endif(COMPILE_TEST_GENERATION_PROGRAMS)

configure_file(${BRAINSFit_SOURCE_DIR}/CTestCustom.ctest ${BRAINSFit_BINARY_DIR}/CTestCustom.ctest COPYONLY)

## Needed to set this for Slicer3 to find the test data directory
#set(BRAINSCommonLib_DATA_DIRS "${PROJECT_SOURCE_DIR}/../BRAINSCommonLi${TestData_DIR}")

set(BRAINSFitTestName BRAINSFitTest_AffineRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --debugLevel 50
)

set(BRAINSFitTestName BRAINSFitTest_AffineRotationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_AffineScaleMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initializeTransformMode useMomentsAlign
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/scale.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/scale.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_AffineScaleNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/scale.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_AffineTranslationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MSE
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initializeTransformMode useMomentsAlign
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/translation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/translation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_AffineTranslationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MSE
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initializeTransformMode useMomentsAlign
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/translation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)


set(BRAINSFitTestName BRAINSFitTest_BSplineAnteScaleRotationRescaleHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500,2500,2500,500
  --numberOfHistogramBins 200
  --splineGridSize 7,5,6
  --numberOfSamples 144000
  --translationScale 250
  --minimumStepLength 0.01,0.003,0.001,0.001
  --outputVolumePixelType short
  --maskProcessingMode ROIAUTO
  --initializeTransformMode useCenterOfHeadAlign
  --transformType Rigid,ScaleVersor3D,Affine,BSpline
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.rescale.rigid.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --maxBSplineDisplacement 7.3
)

set(BRAINSFitTestName BRAINSFitTest_BSplineOnlyRescaleHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 1500
  --numberOfHistogramBins 200
  --splineGridSize 7,5,6
  --numberOfSamples 144000
  --translationScale 250
  --minimumStepLength 0.01
  --outputVolumePixelType short
  --maskProcessingMode ROIAUTO
  --initialTransform DATA{${TestData_DIR}/Initializer_BRAINSFitTest_BSplineAnteScaleRotationRescaleHeadMasks.mat}
  --transformType BSpline
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.rescale.rigid.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --debugLevel 10
  --maxBSplineDisplacement 7.3
  --projectedGradientTolerance 1e-4
  --costFunctionConvergenceFactor 1e+9
)
## NOTE: 7.3 above was computed explicitly through testing.

set(BRAINSFitTestName BRAINSFitTest_BSplineBSplineRescaleHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 1500
  --numberOfHistogramBins 200
  --splineGridSize 7,5,6
  --numberOfSamples 288000
  --translationScale 250
  --minimumStepLength 0.005
  --outputVolumePixelType short
  --maskProcessingMode ROIAUTO
  --initialTransform DATA{${TestData_DIR}/Initializer_BRAINSFitTest_BSplineOnlyRescaleHeadMasks.mat}
  --transformType BSpline
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.rescale.rigid.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --maxBSplineDisplacement 7.3
)

set(BRAINSFitTestName BRAINSFitTest_BSplineScaleRotationRescaleHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500,2500,500
  --numberOfHistogramBins 200
  --splineGridSize 7,5,6
  --numberOfSamples 144000
  --translationScale 250
  --minimumStepLength 0.01,0.003,0.01
  --outputVolumePixelType short
  --maskProcessingMode ROIAUTO
  --initializeTransformMode useCenterOfHeadAlign
  --transformType Rigid,ScaleVersor3D,BSpline
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.rescale.rigid.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --maxBSplineDisplacement 7.3
)


set(BRAINSFitTestName BRAINSFitTest_BSplineScaleRotationHistogramHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1500
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500,2500,500
  --numberOfHistogramBins 50
  --numberOfMatchPoints 10
  --histogramMatch
  --splineGridSize 7,5,6
  --numberOfSamples 144000
  --translationScale 250
  --minimumStepLength 0.01,0.003,0.01
  --outputVolumePixelType short
  --maskProcessingMode ROIAUTO
  --initializeTransformMode useCenterOfHeadAlign
  --transformType Rigid,ScaleVersor3D,BSpline
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.rescale.rigid.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --maxBSplineDisplacement 7.3
)


set(BRAINSFitTestName BRAINSFitTest_RigidAnisotropicMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1500
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --fixedVolume DATA{${TestData_DIR}/ANON0006_20_T1_dbg_splayed.nii.gz}
  --movingVolume DATA{${TestData_DIR}/ANON0006_20_T1_sag_twisted.nii.gz}
  --fixedVolumeTimeIndex 0 --movingVolumeTimeIndex 0
  --minimumStepLength 0.001 --numberOfSamples 100000 --numberOfIterations 1500
  --numberOfHistogramBins 200
  --transformType Rigid --initializeTransformMode useMomentsAlign --translationScale 1000
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  )

set(BRAINSFitTestName BRAINSFitTest_RigidMedianRotationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Rigid
  --initializeTransformMode useMomentsAlign
  --medianFilterSize 1,1,1
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_RigidRotGeomNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 144000
  --translationScale 500
  --minimumStepLength 0.05,0.0050
  --outputVolumePixelType short
  --transformType Rigid,Rigid
  --initializeTransformMode useGeometryAlign
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.geom.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_RigidRotaRotaRotNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.05,0.001,0.00075
  --outputVolumePixelType uchar
  --transformType Rigid,Rigid,Rigid
  --initializeTransformMode useCenterOfHeadAlign
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_RigidRotationHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Rigid
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --maskProcessingMode ROIAUTO
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_RigidRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Rigid
  --initializeTransformMode useMomentsAlign
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_RigidRotationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Rigid
  --interpolationMode WindowedSinc
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)
set(BRAINSFitTestName BRAINSFitTest_RigidRotationNoMasksRiginInPlaceInterp)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Rigid
  --interpolationMode ResampleInPlace
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)


set(BRAINSFitTestName BRAINSFitTest_ScaleRotationRescaleHeadMasksNoInit)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500,2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.01,0.003
  --outputVolumePixelType uchar
  --maskProcessingMode ROIAUTO
  --initializeTransformMode useCenterOfHeadAlign
  --useRigid --useScaleVersor3D
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.rescale.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleSkewVersorRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleSkewVersor3D
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleSkewVersorRotationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleSkewVersor3D
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleSkewVersorScaleMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1200
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleSkewVersor3D
  --initializeTransformMode useMomentsAlign
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/scale.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/scale.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleSkewVersorScaleNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleSkewVersor3D
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/scale.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleTranslationRescaleHeadMasksNoInit)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001,0.0001
  --outputVolumePixelType uchar
  --maskProcessingMode ROIAUTO
  --initializeTransformMode useCenterOfHeadAlign
  --transformType Rigid,ScaleVersor3D
  --permitParameterVariation 0,0,0,1,1,1
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/translation.rescale.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleTranslationRescaleHeadMasksInit)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --maskProcessingMode ROIAUTO
  --initialTransform DATA{${TestData_DIR}/Initializer_BRAINSFitTest_TranslationRescaleHeadMasks.mat}
  --transformType ScaleVersor3D
  --permitParameterVariation 0,0,0,1,1,1,1,1,1
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/translation.rescale.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleVersorRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleVersor3D
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleVersorRotationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleVersor3D
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleVersorScaleMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleVersor3D
  --initializeTransformMode useMomentsAlign
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/scale.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/scale.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_ScaleVersorScaleNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleVersor3D
  --initialTransform DATA{${TestData_DIR}/Initializer_0.05_${BRAINSFitTestName}.mat}
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/scale.test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

if(Test_SignedDistanceData)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSFitTest_SignedDistanceData
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/TEST_output.nii.gz} ${CMAKE_CURRENT_BINARY_DIR}/TEST_output.nii.gz
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode 1
  --fixedVolume DATA{${TestData_DIR}/AtlasBrain_SignedDistance.nii.gz}
  --movingVolume DATA{${TestData_DIR}/Unoriented_RawBrain_SignedDistance.nii.gz}
  --fixedVolumeTimeIndex 0 --movingVolumeTimeIndex 0
  --minimumStepLength 0.001 --numberOfSamples 75000 --numberOfIterations 1500
  --numberOfHistogramBins 200
  --transformType Rigid --initializeTransformMode useMomentsAlign --translationScale 1000
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/TEST.mat
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_output.nii.gz
  )
endif(Test_SignedDistanceData)

set(BRAINSFitTestName BRAINSFitTest_BSplineScaleRotationHeadMasksUShort)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 2505
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MMI
  --backgroundFillValue 10000
  --scaleOutputValues
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500,2500,500
  --splineGridSize 7,5,6
  --numberOfSamples 144000
  --translationScale 250
  --minimumStepLength 0.01,0.003,0.01
  --outputVolumePixelType ushort
  --maskProcessingMode ROIAUTO
  --initializeTransformMode useCenterOfHeadAlign
  --transformType Rigid,ScaleVersor3D,BSpline
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotationUShort.rescale.rigid.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --maxBSplineDisplacement 7.3
)

set(BRAINSFitTestName BRAINSFitTest_BSplineScaleRotationHistogramHeadMasksUShort)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1900
  BRAINSFitTest
  --costMetric MMI
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500,2500,500
  --numberOfHistogramBins 50
  --numberOfMatchPoints 10
  --histogramMatch
  --splineGridSize 7,5,6
  --numberOfSamples 144000
  --translationScale 250
  --minimumStepLength 0.01,0.003,0.01
  --outputVolumePixelType short
  --maskProcessingMode ROIAUTO
  --initializeTransformMode useCenterOfHeadAlign
  --transformType Rigid,ScaleVersor3D,BSpline
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotationUShort.rescale.rigid.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --maxBSplineDisplacement 7.3
)

set(BRAINSFitTestName BRAINSFitTest_NCMetricBrainToItself)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/test.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric NC
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 25
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleVersor3D
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_NCAffineRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 250
  BRAINSFitTest
  --costMetric NC
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --debugLevel 50
)

set(BRAINSFitTestName BRAINSFitTest_NCScaleSkewVersorRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric NC
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleSkewVersor3D
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)


set(BRAINSFitTestName BRAINSFitTest_MSEMetricBrainToItself)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/test.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MSE
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 25
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleVersor3D
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_MSEAffineRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MSE
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --debugLevel 50
)

set(BRAINSFitTestName BRAINSFitTest_MSEScaleSkewVersorRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MSE
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleSkewVersor3D
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_MCMetricBrainToItself)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/test.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MC
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 25
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Rigid --initializeTransformMode useMomentsAlign
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --movingVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_MCMetricMaskBrainToItself)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/testMask.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.testMask.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MC
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 25
  --numberOfHistogramBins 200
  --numberOfSamples 72000
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Rigid --initializeTransformMode useMomentsAlign
  --fixedVolume DATA{${TestData_DIR}/testMask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/testMask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.testMask.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

set(BRAINSFitTestName BRAINSFitTest_MCAffineRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
  BRAINSFitTest
  --costMetric MC
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType Affine
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  --debugLevel 50
)

set(BRAINSFitTestName BRAINSFitTest_MCScaleSkewVersorRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000
  BRAINSFitTest
  --costMetric MC
  --failureExitCode -1 --writeTransformOnFailure
  --numberOfIterations 2500
  --numberOfHistogramBins 200
  --numberOfSamples 131072
  --translationScale 250
  --minimumStepLength 0.001
  --outputVolumePixelType uchar
  --transformType ScaleSkewVersor3D
  --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
  --maskProcessingMode ROI
  --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
  --fixedBinaryVolume DATA{${TestData_DIR}/test_mask.nii.gz}
  --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --movingBinaryVolume DATA{${TestData_DIR}/rotation.test_mask.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
)

# Start of broken metric tests
# Add ", KS (Kappa Statistic for binary images with voxel value=255), MRSD (Mean Reciprocal Square Difference), MIH (Mutual Information Histogram), GD (Gradient Difference), CCH (Correlation Coefficient Histogram), MSEH (Mean Squared Error Histogram), NMIH (Normalized Mutual Information Histogram)" to the end of the description for costMetric in BRAINSFitPrimary.xml then add KS, MRSD, MIH, GD, CCH, MSEH, and NMIH as enumeration elements.
set(RUN_EXTRA_METRIC_TESTS OFF)
if(RUN_EXTRA_METRIC_TESTS)
  # KullbackLeiblerCompareHistogramImageToImageMetric requires training images / infrastructure
  # which is different from the other metrics, making it much harder to incorporate into the
  # BRAINSFIT framework.  This test and the relevant code has been commented out for now.
  #set(BRAINSFitTestName BRAINSFitTest_KLMetricBrainToItself)
  #add_test(NAME ${BRAINSFitTestName} COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
  #  --compare ${BRAINSCommonLib_DATA_DIRS}/test.nii.gz
  #            ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  #  --compareIntensityTolerance 9
  #  --compareRadiusTolerance 0
  #  --compareNumberOfPixelsTolerance 1000
  #  BRAINSFitTest
  #  --costMetric KL
  #  --failureExitCode -1 --writeTransformOnFailure
  #  --numberOfIterations 25
  #  --numberOfHistogramBins 200
  #  --numberOfSamples 72000
  #  --translationScale 250
  #  --minimumStepLength 0.001
  #  --outputVolumePixelType uchar
  #  --transformType ScaleVersor3D
  #  --fixedVolume ${BRAINSCommonLib_DATA_DIRS}/test.nii.gz
  #  --movingVolume ${BRAINSCommonLib_DATA_DIRS}/test.nii.gz
  #  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
  #  --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  #)

  # The Kappa Statistic metric is meant for binary images which is why
  # the fixed and moving images are binarized before the tranfsorm is calculated
  # then the transform is applied to the input image.  This produces somethign
  # close for gray scale images but in the next test using binary images it fails
  # to converge even when given many iterations (which it shouldn't need, since
  # it's registering the same image to itself).
  set(BRAINSFitTestName BRAINSFitTest_KSMetricBrainToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/test.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric KS
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid --initializeTransformMode useMomentsAlign
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )

  set(BRAINSFitTestName BRAINSFitTest_KSMetricBrainMaskToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/testMask.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.testMask.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric KS
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid
    #--initializeTransformMode useMomentsAlign
    --fixedVolume DATA{${TestData_DIR}/testMask.nii.gz}
    --movingVolume DATA{${TestData_DIR}/testMask.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.testMask.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )


  set(BRAINSFitTestName BRAINSFitTest_KSAffineRotationMasks)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 777
    BRAINSFitTest
    --costMetric KS
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 2500
    --numberOfHistogramBins 200
    --numberOfSamples 131072
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Affine
    --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
    --maskProcessingMode ROI
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --fixedBinaryVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
    --movingBinaryVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
    --debugLevel 50
  )

  set(BRAINSFitTestName BRAINSFitTest_KSScaleSkewVersorRotationMasks)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/${BRAINSFitTestName}.result.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 11
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric KS
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 2500
    --numberOfHistogramBins 200
    --numberOfSamples 131072
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType ScaleSkewVersor3D
    --initialTransform DATA{${TestData_DIR}/BRAINSFitTest_Initializer_RigidRotationNoMasks.mat}
    --maskProcessingMode ROI
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --fixedBinaryVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
    --movingBinaryVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )

  # Fails to converge when registering a brain to itself.
  set(BRAINSFitTestName BRAINSFitTest_MRSDMetricBrainToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/test.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric MRSD
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )

  # Segfault when running this test
  set(BRAINSFitTestName BRAINSFitTest_MIHMetricBrainToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/test.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric MIH
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid --initializeTransformMode useMomentsAlign
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )

  # Fails to converge when registering a brain to itself.
  set(BRAINSFitTestName BRAINSFitTest_GDMetricBrainToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/test.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric GD
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid --initializeTransformMode useMomentsAlign
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )

  # Segfault when running this test
  set(BRAINSFitTestName BRAINSFitTest_CCHMetricBrainToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/test.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric CCH
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid --initializeTransformMode useMomentsAlign
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )

  # Segfault when running this test
  set(BRAINSFitTestName BRAINSFitTest_MSEHMetricBrainToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/test.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric MSEH
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid --initializeTransformMode useMomentsAlign
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )

  # Segfault when running this test
  set(BRAINSFitTestName BRAINSFitTest_NMIHMetricBrainToItself)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${BRAINSFitTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSFitTestDriver>
    --compare DATA{${TestData_DIR}/test.nii.gz}
              ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
    BRAINSFitTest
    --costMetric NMIH
    --failureExitCode -1 --writeTransformOnFailure
    --numberOfIterations 25
    --numberOfHistogramBins 200
    --numberOfSamples 72000
    --translationScale 250
    --minimumStepLength 0.001
    --outputVolumePixelType uchar
    --transformType Rigid --initializeTransformMode useMomentsAlign
    --fixedVolume DATA{${TestData_DIR}/test.nii.gz}
    --movingVolume DATA{${TestData_DIR}/test.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.test.nii.gz
    --outputTransform ${CMAKE_CURRENT_BINARY_DIR}/${BRAINSFitTestName}.mat
  )
endif(RUN_EXTRA_METRIC_TESTS)


ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
