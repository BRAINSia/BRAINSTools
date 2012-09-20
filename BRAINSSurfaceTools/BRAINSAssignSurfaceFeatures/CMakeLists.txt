StandardBRAINSBuildMacro(NAME BRAINSAssignSurfaceFeatures TARGET_LIBRARIES ${VTK_LIBRARIES})

if (BUILD_TESTING)

set(BRAINSAssignSurfaceFeaturesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BRAINSAssignSurfaceFeatures)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME TEST_BRAINSAssignSurfaceFeatures
  COMMAND ${LAUNCH_EXE} ${BRAINSAssignSurfaceFeaturesProg}
          --inputSurfaceFile DATA{${TestData_DIR}/0036195_SurfaceGeneration.vtk}
          --outerSurfaceFile DATA{${TestData_DIR}/0036195_leftTissueClass130.vtk}
          --outputSurfaceFile 0036195_SurfaceFeaturesTest.vtk
          --distanceToPC_AP
          --distanceToPC_IS
          --PC -129.5,-178.75,127.97
          --distanceToHull
          --corticalThickness
          --maxThickness 7.0
          --curvature
          --curvatureType Mean )

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME COMPARE_BRAINSAssignSurfaceFeatures
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile 0036195_SurfaceFeaturesTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036195_leftSurfaceFeatures.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_BRAINSAssignSurfaceFeatures PROPERTIES DEPENDS TEST_BRAINSAssignSurfaceFeatures )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

endif ()
