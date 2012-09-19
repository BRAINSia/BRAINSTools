StandardBRAINSBuildMacro(NAME SmoothingMeshScalars TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name SmoothingMeshScalars)

set(SmoothingMeshScalarsProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/SmoothingMeshScalars)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceSmoothingMeshScalars
  COMMAND ${LAUNCH_EXE} ${SmoothingMeshScalarsProg}
          --inputSurfaceFile DATA{${TestData_DIR}/Piecewised_IC4.vtk}
          --outputSurfaceFile SmoothedScalarTest.vtk )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceSmoothingMeshScalars
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile SmoothedScalarTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SmoothedScalar_IC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceSmoothingMeshScalars PROPERTIES DEPENDS TEST_SurfaceSmoothingMeshScalars )

## - ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif ()
