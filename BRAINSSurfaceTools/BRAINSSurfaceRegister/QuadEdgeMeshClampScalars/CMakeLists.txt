StandardBRAINSBuildMacro(NAME QuadEdgeMeshClampScalars TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name QuadEdgeMeshClampScalars)

set(QuadEdgeMeshClampScalarsProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/QuadEdgeMeshClampScalars)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceClampScalars
  COMMAND ${LAUNCH_EXE} ${QuadEdgeMeshClampScalarsProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_AssignArray.vtk}
          --outputSurfaceFile 0036195_ClampScalarsTest.vtk
          --outputMin 10.0
          --outputMax 20.0 )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceClampScalars
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile 0036195_ClampScalarsTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036195_ClampScalars.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceClampScalars PROPERTIES DEPENDS TEST_SurfaceClampScalars )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)
