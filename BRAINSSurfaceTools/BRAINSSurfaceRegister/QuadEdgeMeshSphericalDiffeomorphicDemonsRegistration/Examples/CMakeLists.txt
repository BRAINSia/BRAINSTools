#
# Add here command lines to build examples.
#
SET(QuadEdgeMeshRigidRegistration ${CMAKE_CURRENT_LIST_DIR}/../../QuadEdgeMeshRigidRegistration)
SET(QuadEdgeMeshFieldSmoothingFilters ${CMAKE_CURRENT_LIST_DIR}/../../QuadEdgeMeshFieldSmoothingFilters)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/Testing
  ${QuadEdgeMeshRigidRegistration}/Source/
  ${QuadEdgeMeshFieldSmoothingFilters}/Source/
  )

ADD_EXECUTABLE(QuadEdgeMeshSimilarityCalculator
 QuadEdgeMeshSimilarityCalculator.cxx)
TARGET_LINK_LIBRARIES(QuadEdgeMeshSimilarityCalculator
  ${ITK_LIBRARIES}
)
ADD_EXECUTABLE(ScalarQuadEdgeMeshHistogram ScalarQuadEdgeMeshHistogram.cxx)
TARGET_LINK_LIBRARIES(ScalarQuadEdgeMeshHistogram
  ${ITK_LIBRARIES}
)

IF(USE_VTK)
  ADD_EXECUTABLE(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
   QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1.cxx
   ${QuadEdgeMeshRigidRegistration}/Testing/RegistrationMonitor.cxx)
  TARGET_LINK_LIBRARIES(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
   ITKQuadEdgeMesh ITKNumerics ITKIO
   vtkRendering vtkFiltering vtkIO vtkHybrid)
ELSE(USE_VTK)
  ADD_EXECUTABLE(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
   QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1.cxx)
  TARGET_LINK_LIBRARIES(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
   ITKQuadEdgeMesh ITKNumerics ITKIO)
ENDIF(USE_VTK)

ADD_EXECUTABLE(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter2
 QuadEdgeMeshSphericalDiffeomorphicDemonsFilter2.cxx)
TARGET_LINK_LIBRARIES(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter2
  ${ITK_LIBRARIES})

IF(USE_VTK)
  ADD_EXECUTABLE(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3
   QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3.cxx
   ${QuadEdgeMeshRigidRegistration}/Testing/RegistrationMonitor.cxx)
  TARGET_LINK_LIBRARIES(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3
   ITKQuadEdgeMesh ITKNumerics ITKIO
   vtkRendering vtkFiltering vtkIO vtkHybrid)
ELSE(USE_VTK)
  ADD_EXECUTABLE(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3
   QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3.cxx)
  TARGET_LINK_LIBRARIES(QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3
    ${ITK_LIBRARIES})
ENDIF(USE_VTK)

ADD_EXECUTABLE(ReadFreesurferBinarySurfaceWriteVTKMesh
 ReadFreesurferBinarySurfaceWriteVTKMesh.cxx)
TARGET_LINK_LIBRARIES(ReadFreesurferBinarySurfaceWriteVTKMesh
  ${ITK_LIBRARIES})

ADD_EXECUTABLE(GenerateRegularSphere GenerateRegularSphere.cxx)
TARGET_LINK_LIBRARIES(GenerateRegularSphere
  ${ITK_LIBRARIES})

ADD_EXECUTABLE(ResampleMeshUsingRegularSphere ResampleMeshUsingRegularSphere.cxx)
TARGET_LINK_LIBRARIES(ResampleMeshUsingRegularSphere
  ${ITK_LIBRARIES})

ADD_EXECUTABLE(SmoothingMeshWithScalars SmoothingMeshWithScalars.cxx)
TARGET_LINK_LIBRARIES(SmoothingMeshWithScalars
  ${ITK_LIBRARIES})


IF(USE_VTK)
  ADD_EXECUTABLE(RigidMeshRegistration
    ${QuadEdgeMeshRigidRegistration}/Examples/MeanSquaresMeshToMeshMetricRegistrationWithMonitor1.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/AffineRegistrationMonitor.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/RegistrationMonitor.cxx)
  TARGET_LINK_LIBRARIES(RigidMeshRegistration ITKQuadEdgeMesh ITKIO ITKNumerics
    vtkRendering vtkFiltering vtkIO vtkHybrid)
ELSE(USE_VTK)
  ADD_EXECUTABLE(RigidMeshRegistration
    ${QuadEdgeMeshRigidRegistration}/Examples/MeanSquaresMeshToMeshMetricRegistration1.cxx)
  TARGET_LINK_LIBRARIES(RigidMeshRegistration ITKQuadEdgeMesh ITKIO ITKNumerics)
ENDIF(USE_VTK)

IF(USE_VTK)
  ADD_EXECUTABLE(MultiResolutionRegistration1 MultiResolutionRegistration1.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/AffineRegistrationMonitor.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/RegistrationMonitor.cxx)
  TARGET_LINK_LIBRARIES(MultiResolutionRegistration1 ITKQuadEdgeMesh ITKIO ITKNumerics
    vtkRendering vtkFiltering vtkIO vtkHybrid)
ELSE(USE_VTK)
  ADD_EXECUTABLE(MultiResolutionRegistration1 MultiResolutionRegistration1.cxx)
  TARGET_LINK_LIBRARIES(MultiResolutionRegistration1 ITKQuadEdgeMesh ITKIO ITKNumerics)
ENDIF(USE_VTK)

IF(USE_VTK)
  ADD_EXECUTABLE(MultiResolutionRegistration2 MultiResolutionRegistration2.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/AffineRegistrationMonitor.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/RegistrationMonitor.cxx)
  TARGET_LINK_LIBRARIES(MultiResolutionRegistration2 ITKQuadEdgeMesh ITKIO ITKNumerics
    vtkRendering vtkFiltering vtkIO vtkHybrid)
ELSE(USE_VTK)
  ADD_EXECUTABLE(MultiResolutionRegistration2 MultiResolutionRegistration2.cxx)
  TARGET_LINK_LIBRARIES(MultiResolutionRegistration2 ITKQuadEdgeMesh ITKIO ITKNumerics)
ENDIF(USE_VTK)


IF(USE_VTK)
  ADD_EXECUTABLE(RigidAndDemonsRegistration RigidAndDemonsRegistration.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/AffineRegistrationMonitor.cxx
    ${QuadEdgeMeshRigidRegistration}/Testing/RegistrationMonitor.cxx)
  TARGET_LINK_LIBRARIES(RigidAndDemonsRegistration ITKQuadEdgeMesh ITKIO ITKNumerics
    vtkRendering vtkFiltering vtkIO vtkHybrid)
ELSE(USE_VTK)
  ADD_EXECUTABLE(RigidAndDemonsRegistration RigidAndDemonsRegistration.cxx)
  TARGET_LINK_LIBRARIES(RigidAndDemonsRegistration ITKQuadEdgeMesh ITKIO ITKNumerics)
ENDIF(USE_VTK)

ADD_EXECUTABLE(QuadEdgeMeshSphericalRMSDistanceMetric
 QuadEdgeMeshSphericalRMSDistanceMetric.cxx)
TARGET_LINK_LIBRARIES(QuadEdgeMeshSphericalRMSDistanceMetric
 ITKQuadEdgeMesh ITKNumerics ITKIO)

FUNCTION(TEST_DEPENDENCY TEST1 TEST2)
SET_TESTS_PROPERTIES( ${TEST1} PROPERTIES DEPENDS ${TEST2} )
ENDFUNCTION()


FUNCTION( GENERATE_SPHERE RADIUS SUBDIVISIONS )

ADD_TEST(GenerateIC${SUBDIVISIONS}
  ${CXX_TEST_PATH}/GenerateRegularSphere
  ${TEMP}/IC${SUBDIVISIONS}.vtk
  ${RADIUS}         # Sphere radius
  ${SUBDIVISIONS}   # Number of subdivisions
  )

ENDFUNCTION( GENERATE_SPHERE )


GENERATE_SPHERE( 100.0  0 )
GENERATE_SPHERE( 100.0  1 )
GENERATE_SPHERE( 100.0  2 )
GENERATE_SPHERE( 100.0  3 )
GENERATE_SPHERE( 100.0  4 )
GENERATE_SPHERE( 100.0  5 )
GENERATE_SPHERE( 100.0  6 )
GENERATE_SPHERE( 100.0  7 )



FUNCTION(CONVERT_ANON_DATA ID)

ADD_TEST(anonConvert${ID}a
  ${CXX_TEST_PATH}/ReadFreesurferBinarySurfaceWriteVTKMesh
  ${PATIENT_DATA_ROOT}/sub${ID}/surf/lh.sphere
  ${PATIENT_DATA_ROOT}/sub${ID}/surf/lh.inflated.H
  ${TEMP}/sub${ID}a.vtk
  1 # Normalize the scalar data using robust method
  )

ADD_TEST(anonConvert${ID}b
  ${CXX_TEST_PATH}/ReadFreesurferBinarySurfaceWriteVTKMesh
  ${PATIENT_DATA_ROOT}/sub${ID}/surf/lh.sphere
  ${PATIENT_DATA_ROOT}/sub${ID}/surf/lh.sulc
  ${TEMP}/sub${ID}b.vtk
  1 # Normalize the scalar data using robust method
  )

ADD_TEST(anonConvert${ID}c
  ${CXX_TEST_PATH}/ReadFreesurferBinarySurfaceWriteVTKMesh
  ${PATIENT_DATA_ROOT}/sub${ID}/surf/lh.sphere
  ${PATIENT_DATA_ROOT}/sub${ID}/surf/lh.curv
  ${TEMP}/sub${ID}c.vtk
  1 # Normalize the scalar data using robust method
  )

ENDFUNCTION(CONVERT_ANON_DATA)


FUNCTION( RESAMPLE_MESH INPUT_MESH_ID SUBDIVISIONS)

  ADD_TEST(Resample${INPUT_MESH_ID}_IC${SUBDIVISIONS}
    ${CXX_TEST_PATH}/ResampleMeshUsingRegularSphere
    ${TEMP}/IC${SUBDIVISIONS}.vtk
    ${TEMP}/sub${INPUT_MESH_ID}.vtk
    ${TEMP}/sub${INPUT_MESH_ID}_IC${SUBDIVISIONS}.vtk
    )

ENDFUNCTION( RESAMPLE_MESH )

FUNCTION( SMOOTH_AND_RESAMPLE_MESH INPUT_MESH_ID SUBDIVISIONS1 SUBDIVISIONS2 ITERATIONS )

  ADD_TEST(Smooth${INPUT_MESH_ID}_IC${SUBDIVISIONS1}_${ITERATIONS}
    ${CXX_TEST_PATH}/SmoothingMeshWithScalars
    ${TEMP}/sub${INPUT_MESH_ID}_IC${SUBDIVISIONS1}.vtk
    ${TEMP}/sub${INPUT_MESH_ID}_IC${SUBDIVISIONS1}_SmoothedBy_${ITERATIONS}.vtk
    1.0            # lambda value
    ${ITERATIONS}  # number of iterations
    )

  ADD_TEST(Resample${INPUT_MESH_ID}_IC${SUBDIVISIONS1}_IC${SUBDIVISIONS2}
    ${CXX_TEST_PATH}/ResampleMeshUsingRegularSphere
    ${TEMP}/IC${SUBDIVISIONS2}.vtk
    ${TEMP}/sub${INPUT_MESH_ID}_IC${SUBDIVISIONS1}_SmoothedBy_${ITERATIONS}.vtk
    ${TEMP}/sub${INPUT_MESH_ID}_IC${SUBDIVISIONS2}.vtk
    )

TEST_DEPENDENCY(
  Resample${INPUT_MESH_ID}_IC${SUBDIVISIONS1}_IC${SUBDIVISIONS2}
  Smooth${INPUT_MESH_ID}_IC${SUBDIVISIONS1}_${ITERATIONS}
  )

ENDFUNCTION( SMOOTH_AND_RESAMPLE_MESH )


FUNCTION( RESAMPLE_MESH_ALL_LEVELS INPUT_MESH_ID )

  RESAMPLE_MESH( ${INPUT_MESH_ID}a  4 )  #  IC4: with inflated.H
  RESAMPLE_MESH( ${INPUT_MESH_ID}b  5 )  #  IC5: with sulc
  RESAMPLE_MESH( ${INPUT_MESH_ID}b  6 )  #  IC6: with sulc
  RESAMPLE_MESH( ${INPUT_MESH_ID}c  7 )  #  IC7: with curv

  # Smooth and subsample cascade from IC7 to IC4
  SMOOTH_AND_RESAMPLE_MESH( ${INPUT_MESH_ID}c  7  6  20 )
  SMOOTH_AND_RESAMPLE_MESH( ${INPUT_MESH_ID}c  6  5  10 )
  SMOOTH_AND_RESAMPLE_MESH( ${INPUT_MESH_ID}c  5  4   5 )

ENDFUNCTION( RESAMPLE_MESH_ALL_LEVELS )

FUNCTION( RIGIDLY_REGISTER MESHID1 MESHID2 SUBDIVISIONS )

SET(FIXEDMESH   ${TEMP}/sub${MESHID1}_IC${SUBDIVISIONS}.vtk)
SET(MOVINGMESH  ${TEMP}/sub${MESHID2}_IC${SUBDIVISIONS}.vtk)

ADD_TEST(RigidRegister_${MESHID1}_${MESHID2}_IC${SUBDIVISIONS}
  ${CXX_TEST_PATH}/RigidMeshRegistration
  ${FIXEDMESH}
  ${MOVINGMESH}
  ${TEMP}/sub${MESHID2}_rigid_to_${MESHID1}_IC${SUBDIVISIONS}.vtk
  0 0 1  # Default initial axis
  0      # Default initial angle
  )

ENDFUNCTION( RIGIDLY_REGISTER )

FUNCTION( MULTIRES_REGISTER1 MESHID1 MESHID2 )

SET(FIXEDMESH1   ${TEMP}/sub${MESHID1}a_IC4.vtk)
SET(MOVINGMESH1  ${TEMP}/sub${MESHID2}a_IC4.vtk)

SET(FIXEDMESH2   ${TEMP}/sub${MESHID1}b_IC5.vtk)
SET(MOVINGMESH2  ${TEMP}/sub${MESHID2}b_IC5.vtk)

SET(FIXEDMESH3   ${TEMP}/sub${MESHID1}b_IC6.vtk)
SET(MOVINGMESH3  ${TEMP}/sub${MESHID2}b_IC6.vtk)

SET(FIXEDMESH4   ${TEMP}/sub${MESHID1}c_IC7.vtk)
SET(MOVINGMESH4  ${TEMP}/sub${MESHID2}c_IC7.vtk)

SET(FIXEDMESHFINAL  ${TEMP}/sub${MESHID1}c.vtk)
SET(MOVINGMESHFINAL ${TEMP}/sub${MESHID2}c.vtk)

ADD_TEST(MultiResRegistration1_${MESHID1}_${MESHID2}
  ${CXX_TEST_PATH}/MultiResolutionRegistration1
  ${FIXEDMESH1}
  ${MOVINGMESH1}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_IC4.vtk
  ${FIXEDMESH2}
  ${MOVINGMESH2}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_IC5.vtk
  ${FIXEDMESH3}
  ${MOVINGMESH3}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_IC6.vtk
  ${FIXEDMESH4}
  ${MOVINGMESH4}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_IC7.vtk
  ${FIXEDMESHFINAL}
  ${MOVINGMESHFINAL}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_Final.vtk
  registrationLog_${MESHID2}_${MESHID1}_
  )

ENDFUNCTION( MULTIRES_REGISTER1 )

FUNCTION( MULTIRES_REGISTER2 MESHID1 MESHID2 )

SET(FIXEDMESH1   ${TEMP}/sub${MESHID1}a_IC4.vtk)
SET(MOVINGMESH1  ${TEMP}/sub${MESHID2}a_IC4.vtk)

SET(FIXEDMESH2   ${TEMP}/sub${MESHID1}b_IC5.vtk)
SET(MOVINGMESH2  ${TEMP}/sub${MESHID2}b_IC5.vtk)

SET(FIXEDMESH3   ${TEMP}/sub${MESHID1}b_IC6.vtk)
SET(MOVINGMESH3  ${TEMP}/sub${MESHID2}b_IC6.vtk)

SET(FIXEDMESH4   ${TEMP}/sub${MESHID1}c_IC7.vtk)
SET(MOVINGMESH4  ${TEMP}/sub${MESHID2}c_IC7.vtk)

SET(FIXEDMESHFINAL  ${TEMP}/sub${MESHID1}c.vtk)
SET(MOVINGMESHFINAL ${TEMP}/sub${MESHID2}c.vtk)

ADD_TEST(MultiResRegistration2_${MESHID1}_${MESHID2}
  ${CXX_TEST_PATH}/MultiResolutionRegistration2
  ${FIXEDMESH1}
  ${MOVINGMESH1}
  ${FIXEDMESH2}
  ${MOVINGMESH2}
  ${FIXEDMESH3}
  ${MOVINGMESH3}
  ${FIXEDMESH4}
  ${MOVINGMESH4}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_IC7.vtk
  ${FIXEDMESHFINAL}
  ${MOVINGMESHFINAL}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_Final.vtk
  registrationLog_${MESHID2}_${MESHID1}_
  ${TEMP}/lh_MatlabSphericalDemons_sub${MESHID2}Tosub${MESHID1}.vtk
  4 # number of resolution levels to use
  )

ENDFUNCTION( MULTIRES_REGISTER2 )



FUNCTION( MULTIRES_REGISTER_ONE_FIELD MESHID1 MESHID2 )

SET(FIXEDMESH1   ${TEMP}/sub${MESHID1}c_IC4.vtk)
SET(MOVINGMESH1  ${TEMP}/sub${MESHID2}c_IC4.vtk)

SET(FIXEDMESH2   ${TEMP}/sub${MESHID1}c_IC5.vtk)
SET(MOVINGMESH2  ${TEMP}/sub${MESHID2}c_IC5.vtk)

SET(FIXEDMESH3   ${TEMP}/sub${MESHID1}c_IC6.vtk)
SET(MOVINGMESH3  ${TEMP}/sub${MESHID2}c_IC6.vtk)

SET(FIXEDMESH4   ${TEMP}/sub${MESHID1}c_IC7.vtk)
SET(MOVINGMESH4  ${TEMP}/sub${MESHID2}c_IC7.vtk)

SET(FIXEDMESHFINAL  ${TEMP}/sub${MESHID1}c.vtk)
SET(MOVINGMESHFINAL ${TEMP}/sub${MESHID2}c.vtk)

ADD_TEST(MultiResRegistrationOneField_${MESHID1}_${MESHID2}
  ${CXX_TEST_PATH}/MultiResolutionRegistration2
  ${FIXEDMESH1}
  ${MOVINGMESH1}
  ${FIXEDMESH2}
  ${MOVINGMESH2}
  ${FIXEDMESH3}
  ${MOVINGMESH3}
  ${FIXEDMESH4}
  ${MOVINGMESH4}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_IC7_c.vtk
  ${FIXEDMESHFINAL}
  ${MOVINGMESHFINAL}
  ${TEMP}/sub${MESHID2}_multiresolution_to_${MESHID1}_Final_c.vtk
  registrationLog_${MESHID2}_${MESHID1}_
  ${TEMP}/lh_MatlabSphericalDemons_sub${MESHID2}Tosub${MESHID1}.vtk
  4 # number of resolution levels to use
  )

ENDFUNCTION( MULTIRES_REGISTER_ONE_FIELD )


FUNCTION( RIGID_DEMONS_REGISTER MESHID1 MESHID2 LEVEL RADIUS SIGMAX NUMBEROFITERATIONS)

SET(FIXEDMESH1   ${TEMP}/sub${MESHID1}${LEVEL}.vtk)
SET(MOVINGMESH1  ${TEMP}/sub${MESHID2}${LEVEL}.vtk)

SET(EXPERIMENTID ${MESHID1}_${MESHID2}_${LEVEL})

ADD_TEST(RigidAndDemons_${EXPERIMENTID}
  ${CXX_TEST_PATH}/RigidAndDemonsRegistration
  ${FIXEDMESH1}
  ${MOVINGMESH1}
  ${RADIUS}
  ${SIGMAX}
  ${NUMBEROFITERATIONS}
  ${TEMP}/outputResampledRigid_${EXPERIMENTID}.vtk
  ${TEMP}/outputResampledDemons_${EXPERIMENTID}.vtk
  ${TEMP}/outputDeformedFixedMesh_${EXPERIMENTID}.vtk
  ${TEMP}/DeformationField_${EXPERIMENTID}.vtk
  )

ENDFUNCTION( RIGID_DEMONS_REGISTER )





FUNCTION( DEMONS_REGISTER MESHID1 MESHID2 SUBDIVISIONS )

SET(FIXEDMESH   ${TEMP}/sub${MESHID1}_IC${SUBDIVISIONS}.vtk)
SET(MOVINGMESH  ${TEMP}/sub${MESHID2}_rigid_to_${MESHID1}_IC${SUBDIVISIONS}.vtk)
SET(RESAMPLEDMOVINGMESH  ${TEMP}/sub${MESHID2}_demons_resampled_to_${MESHID1}_IC${SUBDIVISIONS}.vtk)
SET(DEFORMEDFIXEDMESH  ${TEMP}/sub${MESHID1}_demons_deformed_to_${MESHID2}_IC${SUBDIVISIONS}.vtk)

ADD_TEST(DemonsRegister_${MESHID1}_${MESHID2}_IC${SUBDIVISIONS}
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
  ${FIXEDMESH}
  ${MOVINGMESH}
  ${RESAMPLEDMOVINGMESH}
  ${DEFORMEDFIXEDMESH}
  100.0  # Sphere radius
  0.010  # Epsilon equal to 1/(SigmaX^2)
  10.00  # SigmaX  proportional to sqrt(shortest edge). In this mesh short edge is 9.81344
  1.0    # Lambda
   3     # Smoothing Iterations
  15     # maximum number of iterations
  )

ENDFUNCTION( DEMONS_REGISTER )


FUNCTION( COMPARE_MESHES MESHID1 MESHID2 )

SET(REFERENCE_MESH_BASE_NAME
  lh_MatlabSphericalDemons_sub${MESHID1}Tosub${MESHID2})

SET(REFERENCE_MESH
  ${PATIENT_RESULTS_DATA_ROOT}/${REFERENCE_MESH_BASE_NAME}.reg)

ADD_TEST(ConvertResult_${MESHID1}_${MESHID2}
  ${CXX_TEST_PATH}/ReadFreesurferBinarySurfaceWriteVTKMesh
  ${REFERENCE_MESH}
  ${PATIENT_DATA_ROOT}/sub${MESHID2}/surf/lh.curv
  ${TEMP}/${REFERENCE_MESH_BASE_NAME}.vtk
  0  # Do not normalize scalar data
  )


ADD_TEST(CompareRegistration_${MESHID1}_${MESHID2}
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalRMSDistanceMetric
  ${TEMP}/sub${MESHID1}_multiresolution_to_${MESHID2}_Final.vtk
  ${TEMP}/${REFERENCE_MESH_BASE_NAME}.vtk
  ${TEMP}/sub${MESHID1}_multiresolution_to_${MESHID2}_Final_Difference.vtk
  )

ENDFUNCTION( COMPARE_MESHES )


ADD_TEST(DemonsExample1
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC1.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC1.vtk
  ${TEMP}/MeshDemonsFilter1_ResampledMovingMesh_IC1.vtk
  ${TEMP}/MeshDemonsFilter1_DeformedFixedMesh_IC1.vtk
  100.0  # Sphere radius
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
   2     # Smoothing Iterations
  30     # maximum number of iterations
  )

ADD_TEST(DemonsExampleField1
  ${CXX_TEST_PATH}/itkQuadEdgeMeshGenerateDeformationFieldFilterTest1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC1.vtk
  ${TEMP}/MeshDemonsFilter1_DeformedFixedMesh_IC1.vtk
  ${TEMP}/DemonsExampleField1DeformationFieldIC1.vtk
  )
TEST_DEPENDENCY( DemonsExampleField1 DemonsExample1 )


ADD_TEST(DemonsExample2
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC2.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC2.vtk
  ${TEMP}/MeshDemonsFilter1_ResampledMovingMesh_IC2.vtk
  ${TEMP}/MeshDemonsFilter1_DeformedFixedMesh_IC2.vtk
  100.0  # Sphere radius
  0.010  # Epsilon
  10.0   # SigmaX  proportional to sqrt(shortest edge). In this mesh short edge is 39.018
  1.0    # Lambda
   2     # Smoothing Iterations
  30     # maximum number of iterations
  )

ADD_TEST(DemonsExample3
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC3.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC3.vtk
  ${TEMP}/MeshDemonsFilter1_ResampledMovingMesh_IC3.vtk
  ${TEMP}/MeshDemonsFilter1_DeformedFixedMesh_IC3.vtk
  100.0  # Sphere radius
  0.002  # Epsilon equal to 1/(SigmaX^2)
  22.5   # SigmaX  proportional to sqrt(shortest edge). In this mesh short edge is 19.6034
  1.0    # Lambda
   3     # Smoothing Iterations
  50     # maximum number of iterations
  )

ADD_TEST(DemonsExample4
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC4.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC4.vtk
  ${TEMP}/MeshDemonsFilter1_ResampledMovingMesh_IC4.vtk
  ${TEMP}/MeshDemonsFilter1_DeformedFixedMesh_IC4.vtk
  100.0  # Sphere radius
  0.004  # Epsilon equal to 1/(SigmaX^2)
  16.0   # SigmaX  proportional to sqrt(shortest edge). In this mesh short edge is 9.81344
  1.0    # Lambda
   3     # Smoothing Iterations
  50     # maximum number of iterations
  )

ADD_TEST(DemonsExample5
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
  ${TEST_DATA_ROOT}/Input/SimpleSurface1Small.vtk
  ${TEST_DATA_ROOT}/Input/RigidRegistrationSimpleSurface2to1Small.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1Small_ResampledMovingMesh.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1Small_DeformedFixedMesh.vtk
  1.0    # Sphere radius
  0.5    # Epsilon
  0.5    # SigmaX
  1.0    # Lambda
  10     # Smoothing Iterations
  30     # maximum number of iterations
  )

ADD_TEST(DemonsExample6
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter1
  ${TEST_DATA_ROOT}/Input/SimpleSurface1.vtk
  ${TEST_DATA_ROOT}/Input/RigidRegistrationSimpleSurface2to1.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1_ResampledMovingMesh.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1_DeformedFixedMesh.vtk
  100.0  # Sphere radius
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
  10     # Smoothing Iterations
  30     # maximum number of iterations
  )

ADD_TEST(DemonsExample7
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter2
  ${TEST_DATA_ROOT}/Input/SimpleSurface1Small.vtk
  ${TEST_DATA_ROOT}/Input/RigidRegistrationSimpleSurface2to1Small.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1Small_ResampledMovingMesh.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1Small_DeformedFixedMesh.vtk
  1.0    # Sphere radius
  0.5    # Epsilon
  0.5    # SigmaX
  1.0    # Lambda
  10     # Smoothing Iterations
  30     # maximum number of iterations
  )

ADD_TEST(DemonsExample8
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC3.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC3.vtk
  ${TEMP}/MeshDemonsFilter1_ResampledMovingMesh_IC3.vtk
  ${TEMP}/MeshDemonsFilter1_DeformedFixedMesh_IC3.vtk
  100.0  # Sphere radius
  0.002  # Epsilon equal to 1/(SigmaX^2)
  22.5   # SigmaX  proportional to sqrt(shortest edge). In this mesh short edge is 19.6034
  1.0    # Lambda
   3     # Smoothing Iterations
  50     # maximum number of iterations
  3.0    # Metric Significance
  )

ADD_TEST(DemonsExample9
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC4.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC4.vtk
  ${TEMP}/MeshDemonsFilter1_ResampledMovingMesh_IC4.vtk
  ${TEMP}/MeshDemonsFilter1_DeformedFixedMesh_IC4.vtk
  100.0  # Sphere radius
  0.004  # Epsilon equal to 1/(SigmaX^2)
  16.0   # SigmaX  proportional to sqrt(shortest edge). In this mesh short edge is 9.81344
  1.0    # Lambda
   3     # Smoothing Iterations
  50     # maximum number of iterations
  2.0    # Metric Significance
  )

ADD_TEST(DemonsExample10
  ${CXX_TEST_PATH}/QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3
  ${TEST_DATA_ROOT}/Input/SimpleSurface1Small.vtk
  ${TEST_DATA_ROOT}/Input/RigidRegistrationSimpleSurface2to1Small.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1Small_ResampledMovingMesh.vtk
  ${TEMP}/MeshDemonsFilterSimpleSurface2to1Small_DeformedFixedMesh.vtk
  1.0    # Sphere radius
  0.5    # Epsilon
  0.5    # SigmaX
  1.0    # Lambda
  10     # Smoothing Iterations
  30     # maximum number of iterations
  0.01    # Metric Significance
  )

#
# REGISTER( lh 01  15 )
#
# lh_DemonsRegistered_Sub01ToSub15.vtk
#
# lh_moving_to_fixed.vtk
#
# lh_MatlabSphericalDemons_Sub01ToSub15.vtk
# Sub01 is moving
# Sub15 is fixed.
#


IF( PATIENT_DATA_ROOT )

IF( USE_TESTING_SET )
  CONVERT_ANON_DATA(01)
  CONVERT_ANON_DATA(02)
  CONVERT_ANON_DATA(03)
  CONVERT_ANON_DATA(04)
  CONVERT_ANON_DATA(05)
  CONVERT_ANON_DATA(06)
  CONVERT_ANON_DATA(07)
  CONVERT_ANON_DATA(08)
  CONVERT_ANON_DATA(09)
  CONVERT_ANON_DATA(10)

  RESAMPLE_MESH_ALL_LEVELS( 01 )
  RESAMPLE_MESH_ALL_LEVELS( 02 )
  RESAMPLE_MESH_ALL_LEVELS( 03 )
  RESAMPLE_MESH_ALL_LEVELS( 04 )
  RESAMPLE_MESH_ALL_LEVELS( 05 )
  RESAMPLE_MESH_ALL_LEVELS( 06 )
  RESAMPLE_MESH_ALL_LEVELS( 07 )
  RESAMPLE_MESH_ALL_LEVELS( 08 )
  RESAMPLE_MESH_ALL_LEVELS( 09 )
  RESAMPLE_MESH_ALL_LEVELS( 10 )

  RIGIDLY_REGISTER( 01a  02a  4 )
  RIGIDLY_REGISTER( 03a  04a  4 )
  RIGIDLY_REGISTER( 05a  06a  4 )
  RIGIDLY_REGISTER( 07a  08a  4 )
  RIGIDLY_REGISTER( 09a  10a  4 )

  DEMONS_REGISTER( 01a  02a  4 )
  DEMONS_REGISTER( 03a  04a  4 )
  DEMONS_REGISTER( 05a  06a  4 )
  DEMONS_REGISTER( 07a  08a  4 )
  DEMONS_REGISTER( 09a  10a  4 )

# Deprecated
# MULTIRES_REGISTER1( 01 02 )
# MULTIRES_REGISTER1( 03 04 )
# MULTIRES_REGISTER1( 05 06 )
# MULTIRES_REGISTER1( 07 08 )
# MULTIRES_REGISTER1( 09 10 )

  MULTIRES_REGISTER2( 02 01 )
  MULTIRES_REGISTER2( 04 03 )
  MULTIRES_REGISTER2( 06 05 )
  MULTIRES_REGISTER2( 08 07 )
  MULTIRES_REGISTER2( 10 09 )

  COMPARE_MESHES( 01 02 )
  COMPARE_MESHES( 03 04 )
  COMPARE_MESHES( 05 06 )
  COMPARE_MESHES( 07 08 )
  COMPARE_MESHES( 09 10 )

ENDIF( USE_TESTING_SET )

IF( USE_TRAINING_SET )
  CONVERT_ANON_DATA(11)
  CONVERT_ANON_DATA(12)
  CONVERT_ANON_DATA(13)
  CONVERT_ANON_DATA(14)
  CONVERT_ANON_DATA(15)
  CONVERT_ANON_DATA(16)
  CONVERT_ANON_DATA(17)
  CONVERT_ANON_DATA(18)
  CONVERT_ANON_DATA(19)
  CONVERT_ANON_DATA(20)

  RESAMPLE_MESH_ALL_LEVELS( 11 )
  RESAMPLE_MESH_ALL_LEVELS( 12 )
  RESAMPLE_MESH_ALL_LEVELS( 13 )
  RESAMPLE_MESH_ALL_LEVELS( 14 )
  RESAMPLE_MESH_ALL_LEVELS( 15 )
  RESAMPLE_MESH_ALL_LEVELS( 16 )
  RESAMPLE_MESH_ALL_LEVELS( 17 )
  RESAMPLE_MESH_ALL_LEVELS( 18 )
  RESAMPLE_MESH_ALL_LEVELS( 19 )
  RESAMPLE_MESH_ALL_LEVELS( 20 )

  RESAMPLE_MESH( 11b  4 )
  RESAMPLE_MESH( 11b  7 )

  RESAMPLE_MESH( 12b  4 )
  RESAMPLE_MESH( 12b  7 )

  RIGIDLY_REGISTER( 12a  11a  4 )
  RIGIDLY_REGISTER( 14a  13a  4 )
  RIGIDLY_REGISTER( 16a  15a  4 )
  RIGIDLY_REGISTER( 18a  17a  4 )
  RIGIDLY_REGISTER( 20a  19a  4 )

  DEMONS_REGISTER( 12a  11a  4 )
  DEMONS_REGISTER( 14a  13a  4 )
  DEMONS_REGISTER( 16a  15a  4 )
  DEMONS_REGISTER( 18a  17a  4 )
  DEMONS_REGISTER( 20a  19a  4 )

# Deprecated
# MULTIRES_REGISTER1( 12 11 )
# MULTIRES_REGISTER1( 14 13 )
# MULTIRES_REGISTER1( 16 15 )
# MULTIRES_REGISTER1( 18 17 )
# MULTIRES_REGISTER1( 20 19 )

  MULTIRES_REGISTER_ONE_FIELD( 12 11 )

  MULTIRES_REGISTER2( 12 11 )
  MULTIRES_REGISTER2( 14 13 )
  MULTIRES_REGISTER2( 16 15 )
  MULTIRES_REGISTER2( 18 17 )
  MULTIRES_REGISTER2( 20 19 )

  COMPARE_MESHES( 11 12 )
  COMPARE_MESHES( 13 14 )
  COMPARE_MESHES( 15 16 )
  COMPARE_MESHES( 17 18 )
  COMPARE_MESHES( 19 20 )

  RIGID_DEMONS_REGISTER( 12 11 b_IC5  100.0  15.0  30 )
  RIGID_DEMONS_REGISTER( 14 13 b_IC5  100.0  15.0  30 )
  RIGID_DEMONS_REGISTER( 16 15 b_IC5  100.0  15.0  30 )
  RIGID_DEMONS_REGISTER( 18 17 b_IC5  100.0  15.0  30 )
  RIGID_DEMONS_REGISTER( 20 19 b_IC5  100.0  15.0  30 )

  RIGID_DEMONS_REGISTER( 12 11 b_IC6  100.0   9.0  30 )

ENDIF( USE_TRAINING_SET )


ENDIF( PATIENT_DATA_ROOT )

