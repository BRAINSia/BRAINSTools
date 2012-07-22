
include_directories( ${iccdefRegistrationNew_SOURCE_DIR}/Common)


add_executable( ApplyWarpCLPTests  ApplyWarpCLPTests.cxx  )
target_link_libraries(ApplyWarpCLPTests ITKAlgorithms ITKIO)
set(Apply_Warp_TESTS ApplyWarpCLPTests )

add_executable( ICCDEFWarpCLPTests
       iccdefRegistration_New.cxx ICCDEFWarpCLPTests.cxx )
target_link_libraries(ICCDEFWarpCLPTests
                      ${FFTWF_LIB} ITKAlgorithms ITKIO
)

if(USE_DEBUG_IMAGE_VIEWER)
  target_link_libraries(ICCDEFWarpCLPTests ${FFTWF_LIB} ${KWWidgets_LIBRARIES} )
endif(USE_DEBUG_IMAGE_VIEWER)
set(ICCDEF_Warp_TESTS ICCDEFWarpCLPTests )

add_test(ValidateICCDEF_Test1_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${ICCDEF_Warp_TESTS}
  --compare
    ${iccdefRegistrationNew_SOURCE_DIR}/TestData/iccdef_result1.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/test1/forward/deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  ICCDEFWarpTest
    -m ${iccdefRegistrationNew_SOURCE_DIR}/TestData/image1.nii.gz
    -f ${iccdefRegistrationNew_SOURCE_DIR}/TestData/image2.nii.gz
    -o test1
    --useConsistentIntensity
    -v
    -e
    --numberOfHistogramBins 256
    --numberOfMatchPoints 2
    -n 4
    -i 300,150,25,0
    --similarityWeight 5.0
    --outputDisplacementField
    --outputForwardDisplacementFieldVolume  forwardDF.nii.gz
)

add_test(ValidateICCDEF_Test2_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${ICCDEF_Warp_TESTS}
  --compare
    ${iccdefRegistrationNew_SOURCE_DIR}/TestData/iccdef_result2.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/test2/backward/deformedFixedImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  ICCDEFWarpTest
    -m ${iccdefRegistrationNew_SOURCE_DIR}/TestData/image1.nii.gz
    -f ${iccdefRegistrationNew_SOURCE_DIR}/TestData/image2.nii.gz
    -o test2
    --useConsistentIntensity
    -v
    -e
    --numberOfHistogramBins 256
    --numberOfMatchPoints 2
    -n 4
    -i 300,150,25,0
    --medianFilterSize 0.0,0.0,0.0
    --similarityWeight 5.0
)

add_test(ValidateICCDEF_Test3_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${ICCDEF_Warp_TESTS}
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test1/forward/deformedMovingImage.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/test3/forward/deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  ICCDEFWarpTest
    -m ${iccdefRegistrationNew_SOURCE_DIR}/TestData/image1.nii.gz
    -f ${iccdefRegistrationNew_SOURCE_DIR}/TestData/image2.nii.gz
    -o test3
    --useConsistentIntensity
    -v
    -e
    --numberOfHistogramBins 256
    --numberOfMatchPoints 2
    -n 2
    -i 25,0
    --similarityWeight 5.0
    --initialMovingDisplacementFieldVolume  ${iccdefRegistrationNew_BINARY_DIR}/test1/forward/test1_iteration450resolution2500_forward.nii.gz
    --initialFixedDisplacementFieldVolume ${iccdefRegistrationNew_BINARY_DIR}/test1/backward/test1_iteration450resolution2500_backward.nii.gz
)

add_test(ValidateICCDEF_Test4_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${ICCDEF_Warp_TESTS}
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test4/HistogramReferenceFixedImage.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/test4/forward/deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 20
    --compareIntensityTolerance 10
  ICCDEFWarpTest
    -m ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_rotation.input.nii.gz
    -f ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.output.nii.gz
    -o test4
    --useConsistentIntensity
    -v
    -e
    --outputDisplacementField
    --outputJacobianImage
    --outputDisplacement
    --numberOfHistogramBins 256
    --numberOfMatchPoints 2
    -n 1
    -i 1
    --initializeWithTransform ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.mat
    --medianFilterSize 0,0,0
)

add_test(ValidateICCDEF_Test5_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${ICCDEF_Warp_TESTS}
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test5/forward/deformedMovingImage.nii.gz
    ${iccdefRegistrationNew_SOURCE_DIR}/TestData/deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  ICCDEFWarpTest
    -m ${iccdefRegistrationNew_SOURCE_DIR}/TestData/movingImage.nii.gz
    -f ${iccdefRegistrationNew_SOURCE_DIR}/TestData/fixedImage.nii.gz
    -o test5
    --fixedLandmark ${iccdefRegistrationNew_SOURCE_DIR}/TestData/fixedlk.lnd
    --movingLandmark ${iccdefRegistrationNew_SOURCE_DIR}/TestData/movinglk.lnd
    -v
    -e
    --numberOfHistogramBins 256
    --numberOfMatchPoints 2
    -n 3
    -i 16,0,0
    --useConsistentLandmark
    --landmarkWeight -0.04
    --medianFilterSize 0,0,0
)

add_test(ValidateICCDEF_Test6_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${ICCDEF_Warp_TESTS}
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test6/backward/deformedFixedImage.nii.gz
    ${iccdefRegistrationNew_SOURCE_DIR}/TestData/deformedFixedImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  ICCDEFWarpTest
    -m ${iccdefRegistrationNew_SOURCE_DIR}/TestData/movingImage2.nii.gz
    -f ${iccdefRegistrationNew_SOURCE_DIR}/TestData/fixedImage2.nii.gz
    -o test6
    --fixedBinaryVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/fixedImage2_mask.nii.gz
    --movingBinaryVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/movingImage2_mask.nii.gz
    -v
    -e
    --numberOfHistogramBins 256
    --numberOfMatchPoints 2
    -n 2
    -i 5,0
    --useConsistentIntensity
    --initialMovingDisplacementFieldVolume  ${iccdefRegistrationNew_SOURCE_DIR}/TestData/initialMovingDisplacementField.nii.gz
    --initialFixedDisplacementFieldVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/initialFixedDisplacementField.nii.gz
    --similarityWeight 0.7
    --regularizationWeight 50.0
    --inverseWeight 10.0
    --maskProcessingMode BOBF
)


add_test(ValidateApplyWarpTest1_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${Apply_Warp_TESTS}
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test1/forward/deformedMovingImage.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/test4_deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  ApplyWarpTest
    --inputVolume ${iccdefRegistrationNew_BINARY_DIR}/test1/HistogramModifiedMovingImage.nii.gz
    --referenceVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/image2.nii.gz
    --deformationVolume ${iccdefRegistrationNew_BINARY_DIR}/test1/forwardDF.nii.gz
    --orientationRAI --pixelType uchar
    --outputVolume  test4_deformedMovingImage.nii.gz
)

add_test(ValidateApplyWarpTest2_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${Apply_Warp_TESTS}
 --compare
     ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.output.nii.gz
     ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test2.nii.gz
     --compareIntensityTolerance 0
     --compareRadiusTolerance 0
     --compareNumberOfPixelsTolerance 0
   ApplyWarpTest
    --inputVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_rotation.input.nii.gz
    --referenceVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.output.nii.gz
    --outputVolume ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test2.nii.gz
    --pixelType short
    --warpTransform ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.mat
)

add_test(ValidateApplyWarpTest3_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${Apply_Warp_TESTS}
 --compare
     ${iccdefRegistrationNew_SOURCE_DIR}/TestData/ValidateApplyWarpTest3.result.nii.gz
     ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test3.nii.gz
     --compareIntensityTolerance 10
     --compareRadiusTolerance 0
     --compareNumberOfPixelsTolerance 20
   ApplyWarpTest
    --inputVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_rotation.input.nii.gz
    --referenceVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.output.nii.gz
    --outputVolume ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test3.nii.gz
    --pixelType float
    --interpolationMode WindowedSinc
    --warpTransform ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.mat
)

add_test(ValidateApplyWarpTest4_nii ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${Apply_Warp_TESTS}
 --compare
     ${iccdefRegistrationNew_SOURCE_DIR}/TestData/ValidateApplyWarpTest4.result.nii.gz
     ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test4.nii.gz
     --compareIntensityTolerance 0
     --compareRadiusTolerance 0
     --compareNumberOfPixelsTolerance 0
   ApplyWarpTest
    --inputVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_rotation.input.nii.gz
    --referenceVolume ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.output.nii.gz
    --outputVolume ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test4.nii.gz
    --pixelType binary
    --warpTransform ${iccdefRegistrationNew_SOURCE_DIR}/TestData/BRAINSFitTest_RigidRotationMasks.mat
)

set(UNIT_TESTS UnitTests)

set(UnitTest_SRCS
itkVectorFFTWTest.cxx
itkIterativeInverseDisplacementFieldFilterTest.cxx
itkIterativeInverseDisplacementFieldFilterTest2.cxx
)

#add_test(itkVectorFFTWTest ${UNIT_TESTS}
#   --compare ${iccdefRegistrationNew__BINARY_DIR}/transformedImage.nii.gz
#             ${iccdefRegistrationNew__BINARY_DIR}/originalImage.nii.gz
#  --compareNumberOfPixelsTolerance 20
#   itkVectorFFTWTest
#     ${iccdefRegistrationNew__BINARY_DIR}/originalImage.nii.gz
#     ${iccdefRegistrationNew__BINARY_DIR}/transformedImage.nii.gz
#)

#add_executable(UnitTests UnitTests.cxx ${UnitTest_SRCS})
#target_link_libraries(UnitTests
#    ITKAlgorithms ITKIO ${FFTWF_LIB}
#)


#set(UNIT_TESTS ${iccdefRegistrationNew_BINARY_DIR}/bin/UnitTests)
set(UNIT_TESTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/UnitTests)
set(TEMP ${iccdefRegistrationNew_BINARY_DIR})

add_test(itkVectorFFTWTest ${UNIT_TESTS}
   --compare ${TEMP}/transformedImage.nii.gz
             ${TEMP}/originalImage.nii.gz
  --compareNumberOfPixelsTolerance 20
   itkVectorFFTWTest
     ${TEMP}/originalImage.nii.gz
     ${TEMP}/transformedImage.nii.gz
)

add_test(itkIterativeInverseDisplacementFieldFilterTest ${UNIT_TESTS}
   --compare ${TEMP}/multiProcessImage.nii.gz
             ${TEMP}/singleProcessImage.nii.gz
  --compareNumberOfPixelsTolerance 0
  --compareIntensityTolerance 10
   itkIterativeInverseDisplacementFieldFilterTest
     ${TEMP}/singleProcessImage.nii.gz
     ${TEMP}/multiProcessImage.nii.gz
)

add_test(itkIterativeInverseDisplacementFieldFilterTest2 ${UNIT_TESTS}
   --compare ${TEMP}/originalImage.nii.gz
             ${TEMP}/inverse_test2.nii.gz
  --compareNumberOfPixelsTolerance 0
  --compareIntensityTolerance 10
   itkIterativeInverseDisplacementFieldFilterTest2
     ${TEMP}/originalImage.nii.gz
     ${TEMP}/inverse_test2.nii.gz
)

add_executable(UnitTests ${iccdefRegistrationNew_SOURCE_DIR}/unit_test/UnitTests.cxx ${UnitTest_SRCS})
target_link_libraries(UnitTests
    ITKAlgorithms ITKIO ${FFTWF_LIB}
)


