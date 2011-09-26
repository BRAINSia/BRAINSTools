#
# I could check in test data files but I have to write a program to
# build them, so it might as well get run as a percursor test.
add_executable(BRAINSTransformConvertMakeTestFiles
  BRAINSTransformConvertMakeTestFiles.cxx)
target_link_libraries(BRAINSTransformConvertMakeTestFiles BRAINSCommonLib ${ITK_LIBRARIES})

add_executable(ImageCompare ImageCompare.cxx)
target_link_libraries(ImageCompare BRAINSCommonLib ${ITK_LIBRARIES})

set_target_properties(BRAINSTransformConvertMakeTestFiles
  ImageCompare
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

#      <element>VersorRigid</element>
#      <element>ScaleVersor</element>
#      <element>ScaleSkewVersor</element>
#      <element>Affine</element>
#      <element>BSPlineDeformable</element>
#      <element>DisplacementField</element>
#
# Test strategy:
#   foreach transform type
#     convert transform to each higher-dimensional type


#
# program to generate the test files
add_test(NAME BRAINSTransformConvertMakeTestFilesTest
  COMMAND
  ${LAUNCH_EXE} BRAINSTransformConvertMakeTestFiles ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST BRAINSTransformConvertMakeTestFilesTest
  APPEND PROPERTY DEPENDS BRAINSTransformConvertMakeTestFiles)

#
# need to run BRAINSTransformConvert where they were dumped
set(BRAINSTransformConvert "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BRAINSTransformConvert")
set(BRAINSResample "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BRAINSResample")
set(ImageCompare "${CMAKE_CURRENT_BINARY_DIR}/ImageCompare")

#
# use BRAINSTransformConvert to make conversions
set(VersorRigidTargets VersorRigid ScaleVersor ScaleSkewVersor Affine DisplacementField)
set(ScaleVersorTargets ScaleVersor ScaleSkewVersor Affine DisplacementField)
set(ScaleSkewVersorTargets Affine DisplacementField)
set(AffineTargets Affine DisplacementField)
set(BSplineDeformableTargets DisplacementField)

set(TransformTestImage "${CMAKE_CURRENT_BINARY_DIR}/TransformConvertTestImage.nii.gz")

foreach(SourceTransform VersorRigid ScaleVersor ScaleSkewVersor Affine BSplineDeformable)
  foreach(TargetTransform IN LISTS ${SourceTransform}Targets)

    set(TransformConvertTestName Convert${SourceTransform}To${TargetTransform})
    set(SourceTransformName "${CMAKE_CURRENT_BINARY_DIR}/${SourceTransform}Transform.txt")
    set(TargetTransformName "${CMAKE_CURRENT_BINARY_DIR}/${SourceTransform}To${TargetTransform}.txt")
    set(DisplacementVolumeImage "${CMAKE_CURRENT_BINARY_DIR}/${SourceTransform}Displacement.nii.gz")
    #
    # do the conversion
    add_test(NAME ${TransformConvertTestName} COMMAND
      ${LAUNCH_EXE} ${BRAINSTransformConvert}
      --inputTransform ${SourceTransformName}
      --outputTransform ${TargetTransformName}
      --outputTransformType ${TargetTransform}
      --referenceVolume ${TransformTestImage}         # Only needed for DisplacementField transforms
      --displacementVolume ${DisplacementVolumeImage} # Only needed for DisplacementField transforms
      )
    set_property(TEST ${TestName} APPEND PROPERTY DEPENDS BRAINSTransformConvertMakeTestFilesTest)
    set_property(TEST ${TestName} APPEND PROPERTY DEPENDS BRAINSTransformConvert)


    if(NOT "${SourceTransform}" STREQUAL "${TargetTransform}")
      if(NOT "${TargetTransform}" STREQUAL "DisplacementField")
        #
        # Do the backwards conversion which will fail
        set(BadTransformConvertTestName "Bad${TransformConvertTestName}")
        set(BadTransformConvertTargetName ${BadTransformConvertTestName}.txt)
        add_test(NAME ${BadTransformConvertTestName} COMMAND
          ${LAUNCH_EXE} ${BRAINSTransformConvert}
          --inputTransform ${TargetTransformName}
          --outputTransform ${BadTransformConvertTargetName}
          --outputTransformType ${SourceTransform}
          --referenceVolume ${TransformTestImage}
          --displacementVolume ${DisplacementVolumeImage}
          )
        set_property(TEST ${BadTransformConvertTestName}
          APPEND PROPERTY DEPENDS BRAINSTransformConvertMakeTestFilesTest)
        set_property(TEST ${BadTransformConvertTestName}
          APPEND PROPERTY DEPENDS ${TransformConvertTestName})
        set_tests_properties(${BadTransformConvertTestName} PROPERTIES WILL_FAIL ON)
      endif()
    endif(NOT "${SourceTransform}" STREQUAL "${TargetTransform}")

    # resample the test image using BRAINSResample
    # and original transform
    set(ResampleWSourceTestName ResampleWSource${SourceTransform}${TargetTransform})
    set(ResampleWSourceOutputImage ${ResampleWSourceTestName}.nii.gz)

    add_test(NAME ${ResampleWSourceTestName}
      COMMAND ${LAUNCH_EXE} ${BRAINSResample}
      --inputVolume ${TransformTestImage}
      --outputVolume ${ResampleWSourceOutputImage}
      --pixelType short
      --warpTransform ${SourceTransformName})
    set_property(TEST ${ResampleWSourceTestName} APPEND PROPERTY DEPENDS BRAINSResample)
    set_property(TEST ${TestName} APPEND PROPERTY DEPENDS BRAINSTransformConvertMakeTestFilesTest)

    # resample the test image using BRAINSResample
    # and converted transform
    set(ResampleWTargetTestName ResampleWTarget${SourceTransform}${TargetTransform})
    set(ResampleWTargetOutputImage ${ResampleWTargetTestName}.nii.gz)

    # Displacement field takes different command line
    if(NOT "${TargetTransform}" STREQUAL "DisplacementField")
      add_test(NAME ${ResampleWTargetTestName}
        COMMAND ${LAUNCH_EXE} ${BRAINSResample}
        --inputVolume ${TransformTestImage}
        --outputVolume ${ResampleWTargetOutputImage}
        --pixelType short
        --warpTransform ${TargetTransformName})
    else(NOT "${TargetTransform}" STREQUAL "DisplacementField")
      add_test(NAME ${ResampleWTargetTestName}
        COMMAND ${LAUNCH_EXE} ${BRAINSResample}
        --inputVolume ${TransformTestImage}
        --outputVolume ${ResampleWTargetOutputImage}
        --pixelType short
        --deformationVolume ${DisplacementVolumeImage})
    endif(NOT "${TargetTransform}" STREQUAL "DisplacementField")

    #
    # depend on previous test, and BRAINSResample
    set_property(TEST ${ResampleWTargetTestName} APPEND PROPERTY DEPENDS ${TransformConvertTestName})
    set_property(TEST ${ResampleWTargetTestName} APPEND PROPERTY DEPENDS BRAINSResample)

    set(ImageCompareTestName ImageCompare${SourceTransform}${TargetTransform})
    add_test(NAME ${ImageCompareTestName} COMMAND ${LAUNCH_EXE} ${ImageCompare}
      ${ResampleWSourceOutputImage} ${ResampleWTargetOutputImage})

    set_property(TEST ${ImageCompareTestName} APPEND PROPERTY DEPENDS ImageCompare)
    set_property(TEST ${ImageCompareTestName} APPEND PROPERTY DEPENDS ${ResampleWSourceTestName})
    set_property(TEST ${ImageCompareTestName} APPEND PROPERTY DEPENDS ${ResampleWTargetTestName})

  endforeach(TargetTransform)
endforeach(SourceTransform)