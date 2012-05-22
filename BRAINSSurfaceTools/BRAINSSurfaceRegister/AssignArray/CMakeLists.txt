StandardBRAINSBuildMacro(NAME AssignArray TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name AssignArray)

set(AssignArrayProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/AssignArray)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceAssignArray
  COMMAND ${LAUNCH_EXE} ${AssignArrayProg}
          --targetSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --sourceSurfaceFile DATA{${TestData_DIR}/0036195_leftSurfaceFeatures.vtk}
          --outputSurfaceFile 0036195_AssignArrayTest.vtk
          --setArrayName Hull
          --getArrayName distToHull )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceAssignArray
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile 0036195_AssignArrayTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036195_AssignArray.vtk} 
          --vertexLocation 
          --scalarArray 
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceAssignArray PROPERTIES DEPENDS TEST_SurfaceAssignArray )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)