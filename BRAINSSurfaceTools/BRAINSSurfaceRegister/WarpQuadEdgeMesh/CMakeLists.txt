StandardBRAINSBuildMacro(NAME WarpQuadEdgeMesh TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES})

if (BUILD_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceWarpQuadEdgeMesh
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:WarpQuadEdgeMesh>
          --fixedMeshFile DATA{${TestData_DIR}/0036195_left_sphere.vtk}
          --movingMeshFile DATA{${TestData_DIR}/0036896_left_sphereIC6.vtk}
          --deformedMeshFile DATA{${TestData_DIR}/deformed0036195.vtk}
          --outputMeshFile warped0036896HullTest.vtk)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceWarpQuadEdgeMesh
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile warped0036896HullTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/warped0036896Hull.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceWarpQuadEdgeMesh PROPERTIES DEPENDS TEST_SurfaceWarpQuadEdgeMesh )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
