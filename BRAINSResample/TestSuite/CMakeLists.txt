###  This is the test suite for the BRAINSResample program.
###  It is intended to test a wide set of image types and image
###  physical space descriptions to ensure that all ITK
###  images are treated properly.

MakeTestDriverFromSEMTool(BRAINSResample BRAINSResampleTest.cxx)

## Should provide exactly the same result as BRAINSFitTest_AffineRotationMasks
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleTest4_nii
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
  --compare
  DATA{${TestData_DIR}/BRAINSFitTest_AffineRotationMasks.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test4.nii.gz
  --compareIntensityTolerance 13
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1
  BRAINSResampleTest
  --inputVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test4.nii.gz
  --pixelType uchar
  --warpTransform DATA{${TestData_DIR}/BRAINSFitTest_AffineRotationMasks.${XFRM_EXT}}
  )

## Should provide exactly the same result as BRAINSFit BRAINSFitTest_BSplineAnteScaleRotationRescaleHeadMasks test
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleTest5_nii
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
  --compare
  DATA{${TestData_DIR}/BRAINSFitTest_BSplineAnteScaleRotationRescaleHeadMasks.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test5.nii.gz
  --compareIntensityTolerance 13
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 80
  BRAINSResampleTest
  --inputVolume DATA{${TestData_DIR}/rotation.rescale.rigid.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test5.nii.gz
  --pixelType short
  --warpTransform DATA{${TestData_DIR}/BRAINSFitTest_BSplineAnteScaleRotationRescaleHeadMasks.${XFRM_EXT}}
  )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleTest6_nii
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
  --compare
  DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationNoMasks.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test6.nii.gz
  --compareIntensityTolerance 13
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1
  BRAINSResampleTest
  --inputVolume DATA{${TestData_DIR}/rotation.test.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test6.nii.gz
  --pixelType uchar
  --interpolationMode WindowedSinc
  --warpTransform DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationNoMasks.${XFRM_EXT}}
  )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleTest7_nii
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
  --compare
  DATA{${TestData_DIR}/ValidateBRAINSResampleTest7.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test7.nii.gz
  --compareIntensityTolerance 0
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 0
  BRAINSResampleTest
  --inputVolume DATA{${TestData_DIR}/SUBJ_B_small_T1.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/ValidateBRAINSResampleTest7.result.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test7.nii.gz
  --interpolationMode Linear
  --pixelType uchar
  --warpTransform DATA{${TestData_DIR}/applyWarp_BSpline.${XFRM_EXT}}
  )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleTest1_nii
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
  --compare
  DATA{${TestData_DIR}/applyWarp1.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test1.nii.gz
  --compareNumberOfPixelsTolerance 45
  --compareIntensityTolerance 25
  BRAINSResampleTest
  --inputVolume DATA{${TestData_DIR}/SUBJ_B_small_T1.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/applyWarp1.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test1.nii.gz
  --deformationVolume DATA{${TestData_DIR}/OutDefField.nii.gz}
  --pixelType short --interpolationMode NearestNeighbor
  )
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleTest2_nii
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
  --compare
  DATA{${TestData_DIR}/applyWarp2.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test2.nii.gz
  --compareNumberOfPixelsTolerance 20
  --compareIntensityTolerance 25
  BRAINSResampleTest
  --inputVolume DATA{${TestData_DIR}/SUBJ_B_small_T1.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/applyWarp2.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test2.nii.gz
  --deformationVolume DATA{${TestData_DIR}/OutDefField.nii.gz}
  --interpolationMode BSpline --pixelType float
  )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleTest3_nii
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
  --compare
  DATA{${TestData_DIR}/orientedImagesDemons_test.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test3.nii.gz
  --compareNumberOfPixelsTolerance 20
  --compareIntensityTolerance 25
  BRAINSResampleTest
  --inputVolume DATA{${TestData_DIR}/ANON0006_20_T1_sag_twisted.nii.gz}
  --referenceVolume DATA{${TestData_DIR}/orientedImagesDemons_test.nii.gz}
  --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test3.nii.gz
  --deformationVolume DATA{${TestData_DIR}/OutDefField_orientedImage.nii.gz}
  )

if(0) ##This functionality does not work yet
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ValidateBRAINSResampleInvertTest
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSResampleTestDriver>
    --compare
    DATA{${TestData_DIR}/BRAINSResampleTest_Invert.output.nii.gz}
    ${CMAKE_CURRENT_BINARY_DIR}/InvertTransform_Test.nii.gz
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
    BRAINSResampleTest
    --invertTransform
    --inputVolume DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.output.nii.gz}
    --referenceVolume DATA{${TestData_DIR}/BRAINSFitTest_rotation.input.nii.gz}
    --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/applyWarp_test4.nii.gz
    --pixelType short
    --warpTransform DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.${XFRM_EXT}}
    )
endif()

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
