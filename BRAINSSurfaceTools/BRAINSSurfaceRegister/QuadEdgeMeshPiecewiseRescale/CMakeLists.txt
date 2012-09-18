StandardBRAINSBuildMacro(NAME QuadEdgeMeshPiecewiseRescale TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name QuadEdgeMeshPiecewiseRescale)

set(PiecewiseRescaleProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/QuadEdgeMeshPiecewiseRescale)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfacePiecewiseRescale
  COMMAND ${LAUNCH_EXE} ${PiecewiseRescaleProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_sphereIC4.vtk}
          --outputSurfaceFile Piecewised_IC4Test.vtk )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfacePiecewiseRescale
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile Piecewised_IC4Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/Piecewised_IC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfacePiecewiseRescale PROPERTIES DEPENDS TEST_SurfacePiecewiseRescale )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif ()
