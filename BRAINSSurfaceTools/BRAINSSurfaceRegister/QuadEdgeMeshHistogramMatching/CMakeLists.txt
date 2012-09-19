StandardBRAINSBuildMacro(NAME QuadEdgeMeshHistogramMatching  TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name QuadEdgeMeshHistogramMatching)

set(HistogramMatchingProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/QuadEdgeMeshHistogramMatching)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceHistogramMatching
  COMMAND ${LAUNCH_EXE} ${HistogramMatchingProg}
          --refSurfaceFile DATA{${TestData_DIR}/0036195_left_sphereIC4.vtk}
          --inputSurfaceFile DATA{${TestData_DIR}/0036896_left_sphereIC4.vtk}
          --outputSurfaceFile 0036896_to_0036195_hmIC4Test.vtk )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceHistogramMatching
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile 0036896_to_0036195_hmIC4Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036896_to_0036195_hmIC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceHistogramMatching PROPERTIES DEPENDS TEST_SurfaceHistogramMatching )

## - ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif ()
