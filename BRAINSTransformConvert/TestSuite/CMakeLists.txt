#
# I could check in test data files but I have to write a program to
# build them, so it might as well get run as a percursor test.
add_executable(BRAINSTransformConvertMakeTestFiles
  BRAINSTransformConvertMakeTestFiles.cxx)
target_link_libraries(BRAINSTransformConvertMakeTestFiles BRAINSCommonLib ${ITK_LIBRARIES})

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

add_dependencies(BRAINSTransformConvertMakeTestFiles
  BRAINSTransformConvert
  BRAINSResample)

set_target_properties(BRAINSTransformConvertMakeTestFiles
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

add_test(BRAINSTransformConvertMakeTestFiles
  ${LAUNCH_EXE} BRAINSTransformConvertMakeTestFiles ${CMAKE_CURRENT_BINARY_DIR})