StandardBRAINSBuildMacro(NAME SurfaceLabelCleanUp
  ADDITIONAL_SRCS vtkMaskLabel.cxx vtkPolyDataConnectivityIDFilter.cxx
  TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name SurfaceLabelCleanUp)

set(SurfaceLabelCleanUpProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/SurfaceLabelCleanUp)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceLabelCleanUp
  COMMAND ${LAUNCH_EXE} ${SurfaceLabelCleanUpProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_49Labels.vtk}
          --outputSurfaceFile SurfaceLabelCleanUpTest.vtk )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceLabelCleanUp
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile SurfaceLabelCleanUpTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SurfaceLabelCleanUp.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceLabelCleanUp PROPERTIES DEPENDS TEST_SurfaceLabelCleanUp )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)