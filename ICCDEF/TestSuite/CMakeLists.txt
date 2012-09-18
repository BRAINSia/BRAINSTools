include_directories(${${LOCAL_RROJECT_NAME}_SOURCE_DIR})

MakeTestDriverFromSEMTool(ApplyWarp ApplyWarpTest.cxx)
MakeTestDriverFromSEMTool(CreateMask CreateMaskTest.cxx)
MakeTestDriverFromSEMTool(AverageBrainGenerator AverageBrainGeneratorTest.cxx)
MakeTestDriverFromSEMTool(iccdefRegistration_New iccdefRegistration_NewTest.cxx)


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
#    ${ITK_LIBRARIES} ${FFTWF_LIB}
#)


#set(UNIT_TESTS ${iccdefRegistrationNew_BINARY_DIR}/bin/UnitTests)
set(UNIT_TESTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/UnitTests)
set(TEMP ${iccdefRegistrationNew_BINARY_DIR})

add_executable(UnitTests UnitTests.cxx ${UnitTest_SRCS})
target_link_libraries(UnitTests ${ITK_LIBRARIES} ${FFTWF_LIB})


ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateICCDEF_Test1_nii
  COMMAND iccdefRegistration_NewTestDriver
  --compare
    DATA{${TestData_DIR}/iccdef_result1.nii.gz}
    ${iccdefRegistrationNew_BINARY_DIR}/test1/forward/deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
   iccdefRegistration_NewTest
    -m DATA{${TestData_DIR}/image1.nii.gz}
    -f DATA{${TestData_DIR}/image2.nii.gz}
    --outputDirectory ${iccdefRegistrationNew_BINARY_DIR}
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
    --outputForwardDisplacementFieldVolume  ${iccdefRegistrationNew_BINARY_DIR}/test1/forwardDF.nii.gz
)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateApplyWarpTest1_nii
  COMMAND ApplyWarpTestDriver
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test1/forward/deformedMovingImage.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/ValidateApplyWarpTest1_nii_deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  ApplyWarpTest
    --inputVolume ${iccdefRegistrationNew_BINARY_DIR}/test1/HistogramModifiedMovingImage.nii.gz
    --referenceVolume DATA{${TestData_DIR}/image2.nii.gz}
    --deformationVolume ${iccdefRegistrationNew_BINARY_DIR}/test1/forwardDF.nii.gz
    --orientationRAI --pixelType uchar
    --outputVolume  ${iccdefRegistrationNew_BINARY_DIR}/ValidateApplyWarpTest1_nii_deformedMovingImage.nii.gz
)
set_property(TEST ValidateApplyWarpTest1_nii APPEND PROPERTY
  DEPENDS ValidateICCDEF_Test1_nii)


ExternalData_add_test(${PROJECT_NAME}FetchData NAME  ValidateICCDEF_Test2_nii
  COMMAND iccdefRegistration_NewTestDriver
  --compare
    DATA{${TestData_DIR}/iccdef_result2.nii.gz}
    ${iccdefRegistrationNew_BINARY_DIR}/test2/backward/deformedFixedImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  iccdefRegistration_NewTest
    -m DATA{${TestData_DIR}/image1.nii.gz}
    -f DATA{${TestData_DIR}/image2.nii.gz}
    --outputDirectory ${iccdefRegistrationNew_BINARY_DIR}
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

ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateICCDEF_Test3_nii
  COMMAND iccdefRegistration_NewTestDriver
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test1/forward/deformedMovingImage.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/test3/forward/deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  iccdefRegistration_NewTest
    -m DATA{${TestData_DIR}/image1.nii.gz}
    -f DATA{${TestData_DIR}/image2.nii.gz}
    --outputDirectory ${iccdefRegistrationNew_BINARY_DIR}
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

if(0) ## HACK: This test needs to be fixed!
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateICCDEF_Test4_nii
  COMMAND iccdefRegistration_NewTestDriver
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test4/HistogramReferenceFixedImage.nii.gz
    ${iccdefRegistrationNew_BINARY_DIR}/test4/forward/deformedMovingImage.nii.gz
    --compareNumberOfPixelsTolerance 20
    --compareIntensityTolerance 10
  iccdefRegistration_NewTest
    -m DATA{${TestData_DIR}/BRAINSFitTest_rotation.input.nii.gz}
    -f DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.output.nii.gz}
    --outputDirectory ${iccdefRegistrationNew_BINARY_DIR}
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
    --initializeWithTransform DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.mat}
    --medianFilterSize 0,0,0
)
endif() #HACK: This test needs to be fixed!

ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateICCDEF_Test5_nii
  COMMAND iccdefRegistration_NewTestDriver
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test5/forward/deformedMovingImage.nii.gz
    DATA{${TestData_DIR}/deformedMovingImage.nii.gz}
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  iccdefRegistration_NewTest
    -m DATA{${TestData_DIR}/movingImage.nii.gz}
    -f DATA{${TestData_DIR}/fixedImage.nii.gz}
    --outputDirectory ${iccdefRegistrationNew_BINARY_DIR}
    -o test5
    --fixedLandmark DATA{${TestData_DIR}/fixedlk.lnd}
    --movingLandmark DATA{${TestData_DIR}/movinglk.lnd}
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

if(0) #HACK: This test needs to be fixed!
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateICCDEF_Test6_nii
  COMMAND iccdefRegistration_NewTestDriver
  --compare
    ${iccdefRegistrationNew_BINARY_DIR}/test6/backward/deformedFixedImage.nii.gz
    DATA{${TestData_DIR}/deformedFixedImage.nii.gz}
    --compareNumberOfPixelsTolerance 0
    --compareIntensityTolerance 10
  iccdefRegistration_NewTest
    -m DATA{${TestData_DIR}/movingImage2.nii.gz}
    -f DATA{${TestData_DIR}/fixedImage2.nii.gz}
    --outputDirectory ${iccdefRegistrationNew_BINARY_DIR}
    -o test6
    --fixedBinaryVolume DATA{${TestData_DIR}/fixedImage2_mask.nii.gz}
    --movingBinaryVolume DATA{${TestData_DIR}/movingImage2_mask.nii.gz}
    -v
    -e
    --numberOfHistogramBins 256
    --numberOfMatchPoints 2
    -n 2
    -i 5,0
    --useConsistentIntensity
    --initialMovingDisplacementFieldVolume  DATA{${TestData_DIR}/initialMovingDisplacementField.nii.gz}
    --initialFixedDisplacementFieldVolume DATA{${TestData_DIR}/initialFixedDisplacementField.nii.gz}
    --similarityWeight 0.7
    --regularizationWeight 50.0
    --inverseWeight 10.0
    --maskProcessingMode BOBF
)
endif() #HACK: This test needs to be fixed!

if(0) ## These should be re-written to use BRAINSResample instead of ApplyWarp
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateApplyWarpTest2_nii
  COMMAND ApplyWarpTestDriver
  --compare
     DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.output.nii.gz}
     ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test2.nii.gz
     --compareIntensityTolerance 0
     --compareRadiusTolerance 0
     --compareNumberOfPixelsTolerance 0
   ApplyWarpTest
    --inputVolume DATA{${TestData_DIR}/BRAINSFitTest_rotation.input.nii.gz}
    --referenceVolume DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.output.nii.gz}
    --warpTransform DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.mat}
    --outputVolume ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test2.nii.gz
    --pixelType short
)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateApplyWarpTest3_nii
  COMMAND ApplyWarpTestDriver
  --compare
     DATA{${TestData_DIR}/ValidateApplyWarpTest3.result.nii.gz}
     ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test3.nii.gz
     --compareIntensityTolerance 10
     --compareRadiusTolerance 0
     --compareNumberOfPixelsTolerance 20
   ApplyWarpTest
    --inputVolume DATA{${TestData_DIR}/BRAINSFitTest_rotation.input.nii.gz}
    --referenceVolume DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.output.nii.gz}
    --warpTransform DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.mat}
    --outputVolume ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test3.nii.gz
    --pixelType float
    --interpolationMode WindowedSinc
)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME ValidateApplyWarpTest4_nii
  COMMAND ApplyWarpTestDriver
  --compare
     DATA{${TestData_DIR}/ValidateApplyWarpTest4.result.nii.gz}
     ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test4.nii.gz
     --compareIntensityTolerance 0
     --compareRadiusTolerance 0
     --compareNumberOfPixelsTolerance 0
   ApplyWarpTest
    --inputVolume DATA{${TestData_DIR}/BRAINSFitTest_rotation.input.nii.gz}
    --referenceVolume DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.output.nii.gz}
    --warpTransform DATA{${TestData_DIR}/BRAINSFitTest_RigidRotationMasks.mat}
    --outputVolume ${iccdefRegistrationNew_BINARY_DIR}/applyWarp_test4.nii.gz
    --pixelType binary
)
endif() ## HACK

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

if(0)
endif()
ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
