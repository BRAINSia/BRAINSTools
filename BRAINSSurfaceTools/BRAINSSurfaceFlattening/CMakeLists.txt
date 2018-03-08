StandardBRAINSBuildMacro(NAME BRAINSSurfaceFlattening TARGET_LIBRARIES BRAINSCommonLib
  ${BRAINSSurfaceTools_ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)


ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME TEST_BRAINSSurfaceFlattening0
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSSurfaceFlattening>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_smoothed.vtk}
          --outputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --sphereRadius 100.0 )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME COMPARE_BRAINSSurfaceFlattening0
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --refSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --vertexLocation
          --tolerance 1.0e-7 )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME TEST_BRAINSSurfaceFlattening1
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:BRAINSSurfaceFlattening>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_smoothed.vtk}
          --outputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere1.vtk}
          --sphereRadius 100.0
          --seed
          --seedAxis y
          --rotate
          --maxAxis x )

ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME COMPARE_BRAINSSurfaceFlattening1
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere1.vtk}
          --refSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere1.vtk}
          --vertexLocation
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_BRAINSSurfaceFlattening0 PROPERTIES DEPENDS TEST_BRAINSSurfaceFlattening0 )
set_tests_properties( COMPARE_BRAINSSurfaceFlattening1 PROPERTIES DEPENDS TEST_BRAINSSurfaceFlattening1 )

endif ()
