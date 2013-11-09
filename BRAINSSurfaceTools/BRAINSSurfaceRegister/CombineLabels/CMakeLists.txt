StandardBRAINSBuildMacro(NAME CombineLabels TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES})

if (BUILD_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceCombineLabels
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CombineLabels>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_left_49Labels.vtk}
          --outputSurfaceFile CombineLabel0036195Test.vtk
          --removeLabel 49 )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceCombineLabels
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile CombineLabel0036195Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/CombineLabel0036195.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceCombineLabels PROPERTIES DEPENDS TEST_SurfaceCombineLabels )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
