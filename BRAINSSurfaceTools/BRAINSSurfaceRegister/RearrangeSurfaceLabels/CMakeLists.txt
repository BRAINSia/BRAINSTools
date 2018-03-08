StandardBRAINSBuildMacro(NAME RearrangeSurfaceLabels TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES}
  ${BRAINSSurfaceTools_ITK_LIBRARIES})

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME TEST_RearrangeSurfaceLabels
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:RearrangeSurfaceLabels>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_49Labels.vtk}
          --outputSurfaceFile RearrangedLabelsTest.vtk )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME COMPARE_RearrangeSurfaceLabels
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile RearrangedLabelsTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/RearrangedLabels.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_RearrangeSurfaceLabels PROPERTIES DEPENDS TEST_RearrangeSurfaceLabels )

endif ()
