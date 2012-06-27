StandardBRAINSBuildMacro(NAME ResampleQuadEdgeMesh TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name ResampleQuadEdgeMesh)

set(ResampleQuadEdgeMeshProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ResampleQuadEdgeMesh)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceResampleQuadEdgeMesh
  COMMAND ${LAUNCH_EXE} ${ResampleQuadEdgeMeshProg}
          --inputMeshFile DATA{${TestData_DIR}/SimpleSurface1.vtk}
          --refMeshFile DATA{${TestData_DIR}/0036195_left_sphereIC4.vtk}
          --outputMeshFile SimpleSurfaceIC4Test.vtk )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceResampleQuadEdgeMesh
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile SimpleSurfaceIC4Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SimpleSurfaceIC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceResampleQuadEdgeMesh PROPERTIES DEPENDS TEST_SurfaceResampleQuadEdgeMesh )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)