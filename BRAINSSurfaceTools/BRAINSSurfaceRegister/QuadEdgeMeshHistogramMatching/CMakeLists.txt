StandardBRAINSBuildMacro(NAME QuadEdgeMeshHistogramMatching  TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES}
  ${BRAINSSurfaceTools_ITK_LIBRARIES})

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceHistogramMatching
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:QuadEdgeMeshHistogramMatching>
          --refSurfaceFile DATA{${TestData_DIR}/0036195_left_sphereIC4.vtk}
          --inputSurfaceFile DATA{${TestData_DIR}/0036896_left_sphereIC4.vtk}
          --outputSurfaceFile 0036896_to_0036195_hmIC4Test.vtk )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceHistogramMatching
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile 0036896_to_0036195_hmIC4Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036896_to_0036195_hmIC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceHistogramMatching PROPERTIES DEPENDS TEST_SurfaceHistogramMatching )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
