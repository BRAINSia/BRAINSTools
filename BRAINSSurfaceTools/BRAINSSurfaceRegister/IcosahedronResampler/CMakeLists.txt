StandardBRAINSBuildMacro(NAME IcosahedronResampler TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES}
  ${BRAINSSurfaceTools_ITK_LIBRARIES})

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)


ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME TEST_SurfaceIcosahedronResampler1
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:IcosahedronResampler>
          --inputSurfaceFile DATA{${TestData_DIR}/SimpleSurface1.vtk}
          --outputSurfaceFile SimpleSurface1_test_IC4.vtk
          --resolution 4
          --radius 100 )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME TEST_SurfaceIcosahedronResampler2
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:IcosahedronResampler>
          --inputSurfaceFile DATA{${TestData_DIR}/SimpleSurface2.vtk}
          --outputSurfaceFile SimpleSurface2_test_IC4.vtk
          --resolution 4
          --radius 100 )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME COMPARE_SurfaceIcosahedronResampler1
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile SimpleSurface1_test_IC4.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SimpleSurface1_IC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME COMPARE_SurfaceIcosahedronResampler2
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile SimpleSurface2_test_IC4.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SimpleSurface2_IC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceIcosahedronResampler1 PROPERTIES DEPENDS TEST_SurfaceIcosahedronResampler1 )

set_tests_properties( COMPARE_SurfaceIcosahedronResampler2 PROPERTIES DEPENDS TEST_SurfaceIcosahedronResampler2 )

endif ()
