StandardBRAINSBuildMacro(NAME QuadEdgeMeshPiecewiseRescale TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES}
  ${BRAINSSurfaceTools_ITK_LIBRARIES})

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME TEST_SurfacePiecewiseRescale
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:QuadEdgeMeshPiecewiseRescale>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_sphereIC4.vtk}
          --outputSurfaceFile Piecewised_IC4Test.vtk )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME COMPARE_SurfacePiecewiseRescale
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile Piecewised_IC4Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/Piecewised_IC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfacePiecewiseRescale PROPERTIES DEPENDS TEST_SurfacePiecewiseRescale )

endif ()
