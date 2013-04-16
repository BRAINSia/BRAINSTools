StandardBRAINSBuildMacro(NAME ResampleQuadEdgeMesh TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceResampleQuadEdgeMesh
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:ResampleQuadEdgeMesh>
          --inputMeshFile DATA{${TestData_DIR}/SimpleSurface1.vtk}
          --refMeshFile DATA{${TestData_DIR}/0036195_left_sphereIC4.vtk}
          --outputMeshFile SimpleSurfaceIC4Test.vtk )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceResampleQuadEdgeMesh
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile SimpleSurfaceIC4Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SimpleSurfaceIC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceResampleQuadEdgeMesh PROPERTIES DEPENDS TEST_SurfaceResampleQuadEdgeMesh )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
