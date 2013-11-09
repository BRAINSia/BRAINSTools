StandardBRAINSBuildMacro(NAME RemoveTinyLabels
  ADDITIONAL_SRCS vtkMaskLabel.cxx
  TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES})

if (BUILD_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceRemoveTinyLabels
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:RemoveTinyLabels>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_49Labels.vtk}
          --outputSurfaceFile Removed21And37From49LabelsTest.vtk
          --labelList 21,37 )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceRemoveTinyLabels
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile Removed21And37From49LabelsTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/Removed21And37From49Labels.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceRemoveTinyLabels PROPERTIES DEPENDS TEST_SurfaceRemoveTinyLabels )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
