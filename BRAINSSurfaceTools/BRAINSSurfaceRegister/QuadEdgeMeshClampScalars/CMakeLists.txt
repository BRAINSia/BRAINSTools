StandardBRAINSBuildMacro(NAME QuadEdgeMeshClampScalars TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES}
  ${BRAINSSurfaceTools_ITK_LIBRARIES})

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceClampScalars
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:QuadEdgeMeshClampScalars>
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_AssignArray.vtk}
          --outputSurfaceFile 0036195_ClampScalarsTest.vtk
          --outputMin 10.0
          --outputMax 20.0 )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceClampScalars
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile 0036195_ClampScalarsTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036195_ClampScalars.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceClampScalars PROPERTIES DEPENDS TEST_SurfaceClampScalars )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
