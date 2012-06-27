StandardBRAINSBuildMacro(NAME ProbabilityLabels TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name ProbabilityLabels)

set(ProbabilityLabelsProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ProbabilityLabels)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test( ${project_name}FetchData NAME TEST_SurfaceProbabilityLabels
  COMMAND ${LAUNCH_EXE} ${ProbabilityLabelsProg}
          --inputMeshList DATA{${TestData_DIR}/0436396_left_ManualLabels.vtk}
          --inputMeshList DATA{${TestData_DIR}/0436396_left_ManualLabels.vtk}
          --outputMeshFile ProbabilityLabelsTest.vtk
          --probability
          --mostLikely )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceProbabilityLabels
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile ProbabilityLabelsTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/ProbabilityLabels.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceProbabilityLabels PROPERTIES DEPENDS TEST_SurfaceProbabilityLabels )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)