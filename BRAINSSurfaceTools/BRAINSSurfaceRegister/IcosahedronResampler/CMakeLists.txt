StandardBRAINSBuildMacro(NAME IcosahedronResampler TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_TESTING)
set(project_name IcosahedronResampler)

set(IcosahedronResamplerProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/IcosahedronResampler)
set(CompareSurfacesProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/CompareSurfaces)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceIcosahedronResampler1
  COMMAND ${LAUNCH_EXE} ${IcosahedronResamplerProg}
          --inputSurfaceFile DATA{${TestData_DIR}/SimpleSurface1.vtk}
          --outputSurfaceFile SimpleSurface1_test_IC4.vtk
          --resolution 4
          --radius 100 )

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceIcosahedronResampler2
  COMMAND ${LAUNCH_EXE} ${IcosahedronResamplerProg}
          --inputSurfaceFile DATA{${TestData_DIR}/SimpleSurface2.vtk}
          --outputSurfaceFile SimpleSurface2_test_IC4.vtk
          --resolution 4
          --radius 100 )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceIcosahedronResampler1
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile SimpleSurface1_test_IC4.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SimpleSurface1_IC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_SurfaceIcosahedronResampler2
  COMMAND ${LAUNCH_EXE} ${CompareSurfacesProg}
          --inputSurfaceFile SimpleSurface2_test_IC4.vtk
          --refSurfaceFile DATA{${TestData_DIR}/SimpleSurface2_IC4.vtk}
          --vertexLocation
          --scalarArray
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_SurfaceIcosahedronResampler1 PROPERTIES DEPENDS TEST_SurfaceIcosahedronResampler1 )

set_tests_properties( COMPARE_SurfaceIcosahedronResampler2 PROPERTIES DEPENDS TEST_SurfaceIcosahedronResampler2 )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif ()