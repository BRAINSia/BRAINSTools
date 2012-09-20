
StandardBRAINSBuildMacro(NAME BRAINSSurfaceFlattening TARGET_LIBRARIES ${ITK_LIBRARIES})

if (BUILD_TESTING)

set(BRAINSSurfaceFlatteningProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BRAINSSurfaceFlattening)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME TEST_BRAINSSurfaceFlattening0
  COMMAND ${LAUNCH_EXE} ${BRAINSSurfaceFlatteningProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_smoothed.vtk}
          --outputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --sphereRadius 100.0 )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME COMPARE_BRAINSSurfaceFlattening0
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --refSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --vertexLocation
          --tolerance 1.0e-7 )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME TEST_BRAINSSurfaceFlattening1
  COMMAND ${LAUNCH_EXE} ${BRAINSSurfaceFlatteningProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_smoothed.vtk}
          --outputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere1.vtk}
          --sphereRadius 100.0
          --seed
          --seedAxis y
          --rotate
          --maxAxis x )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME COMPARE_BRAINSSurfaceFlattening1
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere1.vtk}
          --refSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere1.vtk}
          --vertexLocation
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_BRAINSSurfaceFlattening0 PROPERTIES DEPENDS TEST_BRAINSSurfaceFlattening0 )
set_tests_properties( COMPARE_BRAINSSurfaceFlattening1 PROPERTIES DEPENDS TEST_BRAINSSurfaceFlattening1 )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
