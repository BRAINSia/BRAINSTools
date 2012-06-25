StandardBRAINSBuildMacro(NAME SurfaceColor TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name SurfaceColor)

set(SurfaceColorProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/SurfaceColor)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceColor
  COMMAND ${LAUNCH_EXE} ${SurfaceColorProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_49Labels.vtk}
          --labelMapFile DATA{${TestData_DIR}/0036195_left_transformedLabelMap.nii.gz}
          --outputSurfaceFile SurfaceColorTest.vtk
          --numOfLabels 49 )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceColor
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile SurfaceColorTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SurfaceColor.vtk} 
          --vertexLocation 
          --scalarArray 
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceColor PROPERTIES DEPENDS TEST_SurfaceColor )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)