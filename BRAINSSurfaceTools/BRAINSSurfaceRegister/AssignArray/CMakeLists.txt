StandardBRAINSBuildMacro(NAME AssignArray TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES} ${BRAINSSurfaceTools_ITK_LIBRARIES})

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceAssignArray
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:AssignArray>
          --targetSurfaceFile DATA{${TestData_DIR}/0036195_left_5000_sphere0.vtk}
          --sourceSurfaceFile DATA{${TestData_DIR}/0036195_leftSurfaceFeatures.vtk}
          --outputSurfaceFile 0036195_AssignArrayTest.vtk
          --setArrayName Hull
          --getArrayName distToHull )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceAssignArray
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile 0036195_AssignArrayTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036195_AssignArray.vtk}
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceAssignArray PROPERTIES DEPENDS TEST_SurfaceAssignArray )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
