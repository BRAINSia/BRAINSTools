StandardBRAINSBuildMacro(NAME AverageScalarsByResampling TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)

set(AverageScalarsByResamplingProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/AverageScalarsByResampling)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceAverageScalarsByResampling
  COMMAND ${LAUNCH_EXE} ${AverageScalarsByResamplingProg}
          --templateSurfaceWithAverageScalars AverageHullOn0036195Test.vtk
          --sphereWithScalarsList DATA{${TestData_DIR}/0036896_left_sphereIC6.vtk}
          --deformedTemplateMeshList DATA{${TestData_DIR}/deformed0036195.vtk}
          --templateSphereFile DATA{${TestData_DIR}/0036195_left_sphere_distToHull.vtk}
          --templateSurfaceFile DATA{${TestData_DIR}/0036195_left_sphere.vtk}
          --numSubs 1 )

ExternalData_add_test(${PROJECT_NAME}FetchData NAME COMPARE_SurfaceAverageScalarsByResampling
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile AverageHullOn0036195Test.vtk
          --refSurfaceFile DATA{${TestData_DIR}/AverageHullOn0036195.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceAverageScalarsByResampling PROPERTIES DEPENDS TEST_SurfaceAverageScalarsByResampling )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
