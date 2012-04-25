#
# Add here command lines to build examples.
#

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR}/../Testing )

ADD_EXECUTABLE(MeanSquaresMeshToMeshMetricRegistration1 MeanSquaresMeshToMeshMetricRegistration1.cxx)
TARGET_LINK_LIBRARIES(MeanSquaresMeshToMeshMetricRegistration1 ITKQuadEdgeMesh ITKNumerics ITKIO)


ADD_EXECUTABLE(ConvertFreeSurferBinarySurfaceToVTKSurface ConvertFreeSurferBinarySurfaceToVTKSurface.cxx)
TARGET_LINK_LIBRARIES(ConvertFreeSurferBinarySurfaceToVTKSurface ITKQuadEdgeMesh ITKNumerics ITKIO)

ADD_TEST(RigidExample1
  ${CXX_TEST_PATH}/MeanSquaresMeshToMeshMetricRegistration1
  ${TEST_DATA_ROOT}/Input/SimpleSurface1Small.vtk
  ${TEST_DATA_ROOT}/Input/SimpleSurface2Small.vtk
  ${TEMP}/RigidRegistrationSimpleSurface2to1Small.vtk
  0 0 1  # Default initial axis
  0      # Default initial angle
  )

ADD_TEST(RigidExample2
  ${CXX_TEST_PATH}/MeanSquaresMeshToMeshMetricRegistration1
  ${TEST_DATA_ROOT}/Input/SimpleSurface1.vtk
  ${TEST_DATA_ROOT}/Input/SimpleSurface2.vtk
  ${TEMP}/RigidRegistrationSimpleSurface2to1.vtk
  0 0 1  # Default initial axis
  0      # Default initial angle
  )


IF(QuadEdgeMeshRigidRegistration_USE_VTK)

  ADD_EXECUTABLE(MeanSquaresMeshToMeshMetricRegistrationWithMonitor1
    MeanSquaresMeshToMeshMetricRegistrationWithMonitor1.cxx
    ${QuadEdgeMeshRigidRegistration_SOURCE_DIR}/Testing/AffineRegistrationMonitor.cxx
    ${QuadEdgeMeshRigidRegistration_SOURCE_DIR}/Testing/RegistrationMonitor.cxx
    )
  TARGET_LINK_LIBRARIES(MeanSquaresMeshToMeshMetricRegistrationWithMonitor1
    ITKQuadEdgeMesh ITKNumerics ITKIO
    vtkRendering vtkFiltering vtkIO vtkHybrid
    )

ADD_TEST(RigidExample3
  ${CXX_TEST_PATH}/MeanSquaresMeshToMeshMetricRegistrationWithMonitor1
  ${TEST_DATA_ROOT}/Input/SimpleSurface1Small.vtk
  ${TEST_DATA_ROOT}/Input/SimpleSurface2Small.vtk
  ${TEMP}/RigidRegistrationSimpleSurface2to1Small2.vtk
  0 0 1  # Default initial axis
  0      # Default initial angle
  )

ADD_TEST(RigidExample4
  ${CXX_TEST_PATH}/MeanSquaresMeshToMeshMetricRegistrationWithMonitor1
  ${TEST_DATA_ROOT}/Input/SimpleSurface1.vtk
  ${TEST_DATA_ROOT}/Input/SimpleSurface2.vtk
  ${TEMP}/RigidRegistrationSimpleSurface2to1_b.vtk
  0 0 1  # Default initial axis
  0      # Default initial angle
  )


ENDIF(QuadEdgeMeshRigidRegistration_USE_VTK)
