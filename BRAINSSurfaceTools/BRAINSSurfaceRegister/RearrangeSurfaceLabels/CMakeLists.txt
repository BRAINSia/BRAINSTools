StandardBRAINSBuildMacro(NAME RearrangeSurfaceLabels TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name RearrangeSurfaceLabels)

set(RearrangeSurfaceLabelsProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/RearrangeSurfaceLabels)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_RearrangeSurfaceLabels
  COMMAND ${LAUNCH_EXE} ${RearrangeSurfaceLabelsProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_49Labels.vtk}
          --outputSurfaceFile RearrangedLabelsTest.vtk )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_RearrangeSurfaceLabels
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile RearrangedLabelsTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/RearrangedLabels.vtk} 
          --vertexLocation 
          --scalarArray 
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_RearrangeSurfaceLabels PROPERTIES DEPENDS TEST_RearrangeSurfaceLabels )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)
