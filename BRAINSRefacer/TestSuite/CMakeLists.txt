# Tests for BRAINSRefacer
# Authors: Jeffrey Obadal, Alexander Leinoff

#
##
### setup test programs/drivers
##
#

MakeTestDriverFromSEMTool(BRAINSRefacer BRAINSRefacerTest.cxx)

#add all simple CLP programs to this list
set(ALL_PROGS_LIST_BRAINSREFACERTESTS
  BRAINSRefacerBrainDataTest
  BRAINSRefacerCheckDeformationAreaLabelMap
  BRAINSRefacerCheckDeformationAreaLandmarks
  BRAINSRefacerCheckPixelType
  )

#For simple CLP programs without TestDriver wrappers
foreach(prog ${ALL_PROGS_LIST_BRAINSREFACERTESTS})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

### setup pixel types to be tested on

# component types enumerations from itkImageIObase - 2016/7/8
# UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE
# string equivalents will be used from ImageIOBase::GetComponentTypeAsString(IOComonentType t)

set(PIXEL_LIST_BRAINSREFACER
    #unsigned_char
    #char
    #unsigned_short
    short
    #unsigned_int
    #int
    #unsigned_long
    #long
    float
    #double
  )

#
##
### Run all the tests against each pixel type
##
#

set(BRR_inPrefix  "${TestData_DIR}/BRAINSRefacer")
set(BRR_outPrefix "${CMAKE_CURRENT_BINARY_DIR}/Output")
file(MAKE_DIRECTORY ${BRR_outPrefix})
# Loop through all pixel types for each test
foreach(PIXELTYPE ${PIXEL_LIST_BRAINSREFACER} )
  #
  ##
  ### Sanity checks on unit tests and testing data
  ##
  #

  # basic test to make sure the input subject image is the pixel type it's supposed to be
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME BRAINSRefacer_${PIXELTYPE}_CheckSubjectImagePixelType
      COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckPixelType>
      --inputImage
      DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
      --pixelType
      ${PIXELTYPE}
      )

  # basic test to make sure the input subject image is the pixel type it's supposed to be
  set(PARENT_TEST BRAINSRefacer_${PIXELTYPE}_CheckSubjectImagePixelType)
  set(RefacerTestName BRAINSRefacer_${PIXELTYPE}_CheckSubjectNewImagePixelType)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
      COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckPixelType>
      --inputImage
      DATA{${BRR_inPrefix}/${PIXELTYPE}_SubjectNewImage.nii.gz}
      --pixelType
      ${PIXELTYPE}
      )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #
  ##
  ### Landmark tests
  ##
  #

  # basic test to make sure landmarks based deformation will run
  set(PARENT_TEST BRAINSRefacer_${PIXELTYPE}_CheckSubjectNewImagePixelType)
  set(RefacerTestName BRAINSRefacerTest_${PIXELTYPE}_RunTest_Landmarks)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    BRAINSRefacerTest
    --inputImage
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --landmarks
    DATA{${BRR_inPrefix}/Landmarks.fcsv}
    --deformedImageName
    ${BRR_outPrefix}/${PIXELTYPE}_LandmarksRefaced_Subject.nii.gz
    --saveTransform
    --bSplineFileName
    ${BRR_outPrefix}/${PIXELTYPE}_Landmarks_bSpline.h5
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  # Test that landmarks refacing output image type matches input image type
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerBrainDataTest_${PIXELTYPE}_Landmarks_CheckOutputPixelType)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
      COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckPixelType>
      --inputImage
      ${BRR_outPrefix}/${PIXELTYPE}_LandmarksRefaced_Subject.nii.gz
      --pixelType
      ${PIXELTYPE}
      )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  # Test to reload bSpline on the original image and compare the results
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerBrainDataTest_${PIXELTYPE}_RunPreviousBspline_Landmarks_OriginalSubject)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    --compare
    ${BRR_outPrefix}/${PIXELTYPE}_LandmarksRefaced_Subject.nii.gz
    ${BRR_outPrefix}/${PIXELTYPE}_ReloadBSplineOriginal_Landmarks_RefacedSubject.nii.gz
    --compareNumberOfPixelsTolerance 50
    --compareIntensityTolerance 0
    BRAINSRefacerTest
    --inputImage
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --deformedImageName
    ${BRR_outPrefix}/${PIXELTYPE}_ReloadBSplineOriginal_Landmarks_RefacedSubject.nii.gz
    --reusebspline
    --previousBSpline
    ${BRR_outPrefix}/${PIXELTYPE}_Landmarks_bSpline.h5
    ## Previous landmarks?
    --landmarks
    DATA{${BRR_inPrefix}/Landmarks.fcsv}
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #Test to reload bSpline and apply it to a new image
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerBrainDataTest_${PIXELTYPE}_RunPreviousBspline_Landmarks_SubjectNewImage)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    BRAINSRefacerTest
    --inputImage
    DATA{${BRR_inPrefix}/${PIXELTYPE}_SubjectNewImage.nii.gz}
    --deformedImageName
    ${BRR_outPrefix}/${PIXELTYPE}_ReloadBSplineOriginal_Landmarks_RefacedSubjectNewImage.nii.gz
    --reusebspline
    --previousBSpline
    ${BRR_outPrefix}/${PIXELTYPE}_Landmarks_bSpline.h5
    ## Previous landmarks?
    --landmarks
    DATA{${BRR_inPrefix}/Landmarks.fcsv}
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #check that no brain data has changed in the landmarks based deformation
  #note this test requires a reference/ground truth brain label map
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerTest_${PIXELTYPE}_CheckBrainData_Landmarks)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
    --inputOriginal
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --inputRefaced
    ${BRR_outPrefix}/${PIXELTYPE}_LandmarksRefaced_Subject.nii.gz
    --labelmap
    DATA{${BRR_inPrefix}/${PIXELTYPE}_BrainLabels.nii.gz}
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  # check to see that non brain data changed for landmarks
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerTest_${PIXELTYPE}_CheckDeformationArea_Landmarks)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLandmarks>
    --inputOriginal
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --inputRefaced
    ${BRR_outPrefix}/${PIXELTYPE}_LandmarksRefaced_Subject.nii.gz
    --landmarks
    DATA{${BRR_inPrefix}/Landmarks.fcsv}
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_Landmarks)
  set(RefacerTestName BRAINSRefacerTest_${PIXELTYPE}_CheckNonDeformationArea_Landmarks)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLandmarks>
    --inputOriginal
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --inputRefaced
    ${BRR_outPrefix}/${PIXELTYPE}_LandmarksRefaced_Subject.nii.gz
    --landmarks
    DATA{${BRR_inPrefix}/Landmarks.fcsv}
    --checkNonDeformedArea
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #
  ##
  ### labelmap tests
  ##
  #

  # basic test to make sure brain label map based deformation will run
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME BRAINSRefacerTest_${PIXELTYPE}_RunTest_LabelMap
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    BRAINSRefacerTest
    --inputImage
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --labelmap
    DATA{${BRR_inPrefix}/${PIXELTYPE}_BrainLabels.nii.gz}
    --deformedImageName
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_RefacedSubject.nii.gz
    --labelmapMask
    --saveTransform
    --bSplineFileName
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_bSpline.h5
    )

  # Test that landmarks refacing output image type matches input image type
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_LabelMap)
  set(RefacerTestName BRAINSRefacerBrainDataTest_${PIXELTYPE}_LabelMap_CheckOutputPixelType)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
      COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckPixelType>
      --inputImage
      ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_RefacedSubject.nii.gz
      --pixelType
      ${PIXELTYPE}
      )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #Test to reload bSpline and apply it to a new image
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_LabelMap)
  set(RefacerTestName BRAINSRefacerBrainDataTest_${PIXELTYPE}_RunPreviousBspline_LabelMap_NewSubjectImage)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    BRAINSRefacerTest
    --inputImage
    DATA{${BRR_inPrefix}/${PIXELTYPE}_SubjectNewImage.nii.gz}
    --labelmap
    DATA{${BRR_inPrefix}/${PIXELTYPE}_BrainLabels.nii.gz}
    --deformedImageName
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_RefacedNSubjectNewImage.nii.gz
    --labelmapMask
    --reusebspline
    --previousBSpline
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_bSpline.h5
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #check that no brain data has been changed in label map based deformation
  # note this test requires a reference/ground truth brain label map
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_LabelMap)
  set(RefacerTestName BRAINSRefacerTest_${PIXELTYPE}_CheckBrainData_LabelMap)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerBrainDataTest>
    --inputOriginal
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --inputRefaced
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_RefacedSubject.nii.gz
    --labelmap
    DATA{${BRR_inPrefix}/${PIXELTYPE}_BrainLabels.nii.gz}
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  #Test to reload bSpline and apply it to the same image, and compare to see if the output is the same
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_LabelMap)
  set(RefacerTestName BRAINSRefacerBrainDataTest_${PIXELTYPE}_RunPreviousBspline_LabelMap_OriginalSubject)
  ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerTestDriver>
    --compare
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_RefacedSubject.nii.gz
    ${BRR_outPrefix}/${PIXELTYPE}_ReloadBSplineOriginal_LabelMap_RefacedSubject.nii.gz
    --compareNumberOfPixelsTolerance 50
    --compareIntensityTolerance 0
    BRAINSRefacerTest
    --inputImage
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --deformedImageName
    ${BRR_outPrefix}/${PIXELTYPE}_ReloadBSplineOriginal_LabelMap_RefacedSubject.nii.gz
    --reusebspline
    --previousBSpline
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_bSpline.h5
    --labelmapMask
    ## Previous label map
    --labelmap
    DATA{${BRR_inPrefix}/${PIXELTYPE}_BrainLabels.nii.gz}
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  # check to see that non brain data has changed for label map
  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_LabelMap)
  set(RefacerTestName BRAINSRefacerTest_${PIXELTYPE}_CheckDeformationArea_LabelMmap)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLabelMap>
    --inputOriginal
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --inputRefaced
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_RefacedSubject.nii.gz
    --labelmap
    DATA{${BRR_inPrefix}/${PIXELTYPE}_BrainLabels.nii.gz}
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

  set(PARENT_TEST BRAINSRefacerTest_${PIXELTYPE}_RunTest_LabelMap)
  set(RefacerTestName BRAINSRefacerTest_${PIXELTYPE}_CheckNonDeformationArea_LabelMap)
  ExternalData_add_test(${PROJECT_NAME}FetchData NAME ${RefacerTestName}
    COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSRefacerCheckDeformationAreaLabelMap>
    --inputOriginal
    DATA{${BRR_inPrefix}/${PIXELTYPE}_Subject.nii.gz}
    --inputRefaced
    ${BRR_outPrefix}/${PIXELTYPE}_LabelMap_RefacedSubject.nii.gz
    --labelmap
    DATA{${BRR_inPrefix}/${PIXELTYPE}_BrainLabels.nii.gz}
    --checkNonDeformedArea
    )
  set_property(TEST ${RefacerTestName} APPEND PROPERTY DEPENDS ${PARENT_TEST})

# #TODO: Test that the output pixel type is the same as the input

# #TODO: Test that the deformation field is symmetrical

endforeach()
