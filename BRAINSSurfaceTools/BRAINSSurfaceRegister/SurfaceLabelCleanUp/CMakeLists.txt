StandardBRAINSBuildMacro(NAME SurfaceLabelCleanUp
  ADDITIONAL_SRCS vtkMaskLabel.cxx vtkPolyDataConnectivityIDFilter.cxx
  TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES})

if (BUILD_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceLabelCleanUp
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:SurfaceLabelCleanUp>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_49Labels.vtk}
          --outputSurfaceFile SurfaceLabelCleanUpTest.vtk )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceLabelCleanUp
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile SurfaceLabelCleanUpTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SurfaceLabelCleanUp.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceLabelCleanUp PROPERTIES DEPENDS TEST_SurfaceLabelCleanUp )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
