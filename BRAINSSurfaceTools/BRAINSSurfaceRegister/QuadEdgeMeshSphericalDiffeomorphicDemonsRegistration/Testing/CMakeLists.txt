SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
SET(TEST_DATA_ROOT ${CMAKE_SOURCE_DIR}/Data)
SET(TEMP ${CMAKE_BINARY_DIR}/Testing/Temporary)

SET(QuadEdgeMeshRigidRegistration ${PROJECT_SOURCE_DIR}/../QuadEdgeMeshRigidRegistration)
SET(QuadEdgeMeshFieldSmoothingFilters ${PROJECT_SOURCE_DIR}/../QuadEdgeMeshFieldSmoothingFilters)

FIND_PATH(FREESURFER_DIR
  SetUpFreeSurfer.sh
  PATHS
  $ENV{FREESURFER_HOME}
  DOC "Directory where FreeSurfer is installed"
  )

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/Testing
  ${QuadEdgeMeshFieldSmoothingFilters}/Source
  ${QuadEdgeMeshFieldSmoothingFilters}/Testing
  )

#######################
#######################
#
#  Define Executables
#
#######################
#######################


ADD_EXECUTABLE(itkMeshGeneratorHelperTest1 
  ${QuadEdgeMeshFieldSmoothingFilters}/Testing/itkMeshGeneratorHelperTest1.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest2
  ${QuadEdgeMeshFieldSmoothingFilters}/Testing/itkMeshGeneratorHelperTest2.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest3
  ${QuadEdgeMeshFieldSmoothingFilters}/Testing/itkMeshGeneratorHelperTest3.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest3 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest4
  ${QuadEdgeMeshFieldSmoothingFilters}/Testing/itkMeshGeneratorHelperTest4.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest4 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest5 itkMeshGeneratorHelperTest5.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest5 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest6 itkMeshGeneratorHelperTest6.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest6 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkLinearInterpolateDeformationFieldTest1 itkLinearInterpolateDeformationFieldMeshFunctionTest1.cxx )
TARGET_LINK_LIBRARIES(itkLinearInterpolateDeformationFieldTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshSphericalDemonsFilterTest1 itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilterTest1.cxx )
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshSphericalDemonsFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshSphericalDemonsFilterTest2 itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilterTest2.cxx )
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshSphericalDemonsFilterTest2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshSimilarityCalculatorTest itkQuadEdgeMeshSimilarityCalculatorTest.cxx )
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshSimilarityCalculatorTest ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkScalarQuadEdgeMeshToListAdaptorTest1 itkScalarQuadEdgeMeshToListAdaptorTest1.cxx )
TARGET_LINK_LIBRARIES(itkScalarQuadEdgeMeshToListAdaptorTest1 ITKQuadEdgeMesh ITKIO ITKNumerics ITKStatistics)

ADD_EXECUTABLE(itkResampleDestinationPointsQuadEdgeMeshFilterTest1 itkResampleDestinationPointsQuadEdgeMeshFilterTest1.cxx )
TARGET_LINK_LIBRARIES(itkResampleDestinationPointsQuadEdgeMeshFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkDeformQuadEdgeMeshFilterTest1 itkDeformQuadEdgeMeshFilterTest1.cxx )
TARGET_LINK_LIBRARIES(itkDeformQuadEdgeMeshFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkNormalizeScalarsQuadEdgeMeshFilterTest1 itkNormalizeScalarsQuadEdgeMeshFilterTest1.cxx )
TARGET_LINK_LIBRARIES(itkNormalizeScalarsQuadEdgeMeshFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshGenerateDeformationFieldFilterTest1 itkQuadEdgeMeshGenerateDeformationFieldFilterTest1.cxx )
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshGenerateDeformationFieldFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkVersorOperationsTest1 itkVersorOperationsTest1.cxx )
TARGET_LINK_LIBRARIES(itkVersorOperationsTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkNodeScalarGradientCalculatorTest1 itkNodeScalarGradientCalculatorTest1.cxx )
TARGET_LINK_LIBRARIES(itkNodeScalarGradientCalculatorTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkNodeScalarGradientCalculatorTest2 itkNodeScalarGradientCalculatorTest2.cxx )
TARGET_LINK_LIBRARIES(itkNodeScalarGradientCalculatorTest2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkNodeVectorJacobianCalculatorTest1 itkNodeVectorJacobianCalculatorTest1.cxx )
TARGET_LINK_LIBRARIES(itkNodeVectorJacobianCalculatorTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkNodeVectorJacobianCalculatorTest2 itkNodeVectorJacobianCalculatorTest2.cxx )
TARGET_LINK_LIBRARIES(itkNodeVectorJacobianCalculatorTest2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkAssignScalarValuesQuadEdgeMeshFilterTest1 itkAssignScalarValuesQuadEdgeMeshFilterTest1.cxx )
TARGET_LINK_LIBRARIES(itkAssignScalarValuesQuadEdgeMeshFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkTriangleBasisSystemTest1 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkTriangleBasisSystemTest1.cxx )
TARGET_LINK_LIBRARIES(itkTriangleBasisSystemTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkTriangleListBasisSystemCalculatorTest1 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkTriangleListBasisSystemCalculatorTest1.cxx )
TARGET_LINK_LIBRARIES(itkTriangleListBasisSystemCalculatorTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkTriangleBasisSystemCalculatorTest1 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkTriangleBasisSystemCalculatorTest1.cxx )
TARGET_LINK_LIBRARIES(itkTriangleBasisSystemCalculatorTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkLinearInterpolateMeshFunctionTest2 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkLinearInterpolateMeshFunctionTest2.cxx )
TARGET_LINK_LIBRARIES(itkLinearInterpolateMeshFunctionTest2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkLinearInterpolateMeshFunctionTest1 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkLinearInterpolateMeshFunctionTest1.cxx )
TARGET_LINK_LIBRARIES(itkLinearInterpolateMeshFunctionTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkFreeSurferBinarySurfaceReaderTest1 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkFreeSurferBinarySurfaceReaderTest1.cxx )
TARGET_LINK_LIBRARIES(itkFreeSurferBinarySurfaceReaderTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkFreeSurferBinarySurfaceReaderTest2 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkFreeSurferBinarySurfaceReaderTest2.cxx )
TARGET_LINK_LIBRARIES(itkFreeSurferBinarySurfaceReaderTest2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkPointLocator2Test1 
  ${QuadEdgeMeshRigidRegistration}/Testing/itkPointLocator2Test1.cxx )
TARGET_LINK_LIBRARIES(itkPointLocator2Test1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkWarpQuadEdgeMeshFilterTest itkWarpQuadEdgeMeshFilterTest.cxx )
TARGET_LINK_LIBRARIES(itkWarpQuadEdgeMeshFilterTest ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkRescaleScalarsQuadEdgeMeshFilterTest itkRescaleScalarsQuadEdgeMeshFilterTest.cxx)
TARGET_LINK_LIBRARIES(itkRescaleScalarsQuadEdgeMeshFilterTest ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkHistogramMatchingQuadEdgeMeshFilterTest itkHistogramMatchingQuadEdgeMeshFilterTest.cxx)
TARGET_LINK_LIBRARIES(itkHistogramMatchingQuadEdgeMeshFilterTest ITKQuadEdgeMesh ITKIO ITKNumerics ITKStatistics)

IF(USE_VTK)

  ADD_EXECUTABLE(ScalarMeshRender  ScalarMeshRender.cxx)
  TARGET_LINK_LIBRARIES(ScalarMeshRender 
    vtkRendering vtkFiltering vtkIO vtkHybrid
    ) 
    
  ADD_EXECUTABLE(MeshRender  MeshRender.cxx)
  TARGET_LINK_LIBRARIES(MeshRender vtkRendering vtkFiltering vtkIO vtkHybrid) 
    
  ADD_EXECUTABLE(VectorMeshRender  VectorMeshRender.cxx)
  TARGET_LINK_LIBRARIES(VectorMeshRender 
    vtkRendering vtkFiltering vtkIO vtkHybrid
    ) 

ENDIF(USE_VTK)



###################
###################
#
#  Define Tests
#
###################
###################

ADD_TEST(MeshGeneratorHelperTest1
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest1
  ${TEMP}/fixedMesh01.vtk
  ${TEMP}/movingMesh01.vtk
  )

ADD_TEST(MeshGeneratorHelperTest2
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest2
  ${TEMP}/meshWithVectors.vtk
  )

ADD_TEST(MeshGeneratorHelperTest3
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest3
  ${TEMP}/meshWithVectorsIC1.vtk
  1      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(MeshGeneratorHelperTest4
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest3
  ${TEMP}/meshWithVectorsIC4.vtk
  4      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(MeshGeneratorHelperTest5
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest3
  ${TEMP}/meshWithVectorsIC7.vtk
  7      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(MeshGeneratorHelperTest6
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest4
  ${TEMP}/meshWithScalarsIC1.vtk
  1      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(MeshGeneratorHelperTest7
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest4
  ${TEMP}/meshWithScalarsIC4.vtk
  4      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(MeshGeneratorHelperTest8
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest5
  ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk
  4      # resolution
  100.0  # Scale (radius)
  3.0    # Phi multiplier
  )

ADD_TEST(MeshGeneratorHelperTest9
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest5
  ${TEMP}/movingMeshWithScalarsGaussianIC4.vtk
  4      # resolution
  100.0  # Scale (radius)
  2.0    # Phi multiplier
  )
  
ADD_TEST(MeshGeneratorHelperTest10
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest6
  ${TEMP}/meshWithoutScalarsGaussianIC4.vtk
  4      # resolution
  100.0  # Scale (radius)
  )
  
ADD_TEST(LinearInterpolateDeformationFieldMeshFunctionTest1
  ${CXX_TEST_PATH}/itkLinearInterpolateDeformationFieldTest1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC1.vtk
  )
SET_TESTS_PROPERTIES( LinearInterpolateDeformationFieldMeshFunctionTest1 PROPERTIES DEPENDS MeshGeneratorHelperTest2)

ADD_TEST(DemonsFilterTest1
  ${CXX_TEST_PATH}/itkQuadEdgeMeshSphericalDemonsFilterTest1
  ${TEMP}/meshWithScalarsIC1.vtk
  ${TEMP}/meshWithScalarsIC1.vtk
  ${TEMP}/meshWithScalarsIC1RegisteredToItself.vtk
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
  1      # Smoothing Iterations
  1      # maximum number of iterations
  )

ADD_TEST(DemonsFilterTest2
  ${CXX_TEST_PATH}/itkQuadEdgeMeshSphericalDemonsFilterTest1
  ${TEMP}/meshWithScalarsIC4.vtk
  ${TEMP}/meshWithScalarsIC4.vtk
  ${TEMP}/meshWithScalarsIC4RegisteredToItself.vtk
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
  1      # Smoothing Iterations
  1      # maximum number of iterations
  )

ADD_TEST(DemonsFilterTest3
  ${CXX_TEST_PATH}/itkQuadEdgeMeshSphericalDemonsFilterTest1
  ${TEMP}/fixedMesh01.vtk
  ${TEMP}/movingMesh01.vtk
  ${TEMP}/movingMesh01RegisteredToFixedMesh01.vtk
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
  1      # Smoothing Iterations
  1      # maximum number of iterations
  )

ADD_TEST(DemonsFilterTest4
  ${CXX_TEST_PATH}/itkQuadEdgeMeshSphericalDemonsFilterTest1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC1.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC1.vtk
  ${TEMP}/movingMeshWithScalarsIC1registeredToFixedMeshWithScalarsIC1.vtk
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
  10     # Smoothing Iterations
  15     # maximum number of iterations
  )

ADD_TEST(DemonsFilterTest5
  ${CXX_TEST_PATH}/itkQuadEdgeMeshSphericalDemonsFilterTest1
  ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/movingMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/registeredMeshWithScalarsGaussianIC4.vtk
  0.01     # Epsilon
  100.0    # SigmaX
  1.0      # Lambda
  2        # Smoothing Iterations
  20      # maximum number of iterations
  )
  
ADD_TEST(DemonsFilterTest6
  ${CXX_TEST_PATH}/itkQuadEdgeMeshSphericalDemonsFilterTest2
  ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/movingMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/deformedMeshWithScalarsGaussianIC4.vtk
  0.01     # Epsilon
  100.0    # SigmaX
  1.0      # Lambda
  2        # Smoothing Iterations
  20      # maximum number of iterations
  )

ADD_TEST(ResampleDestinationFilterTest1
  ${CXX_TEST_PATH}/itkResampleDestinationPointsQuadEdgeMeshFilterTest1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC3.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC3.vtk
  ${TEMP}/deformedFixedMeshWithScalarsIC3InResamplingTest.vtk
  ${TEMP}/meshWithScalarsIC4.vtk
  ${TEMP}/deformedFixedMeshWithScalarsIC3ResampledToIC4.vtk
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
  1      # Smoothing Iterations
  5      # maximum number of iterations
  )
  
ADD_TEST(DeformFilterTest1
  ${CXX_TEST_PATH}/itkDeformQuadEdgeMeshFilterTest1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC3.vtk
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC3.vtk
  ${TEMP}/deformedFixedMeshWithScalarsIC3.vtk
  ${TEMP}/meshWithScalarsIC4.vtk
  ${TEMP}/deformedFixedMeshWithScalarsIC4.vtk
  1.0    # Epsilon
  1.0    # SigmaX
  1.0    # Lambda
  1      # Smoothing Iterations
  5      # maximum number of iterations
  )
  
ADD_TEST(ComputeDeformationFieldMeshTest1
  ${CXX_TEST_PATH}/itkQuadEdgeMeshGenerateDeformationFieldFilterTest1
  ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/deformedMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/fixedMeshWithDeformationFieldIC4.vtk
  )
SET_TESTS_PROPERTIES( ComputeDeformationFieldMeshTest1 PROPERTIES DEPENDS DemonsFilterTest6)

ADD_TEST(VersorOperationsTest1
  ${CXX_TEST_PATH}/itkVersorOperationsTest1
  )

ADD_TEST(NodeScalarGradientCalculatorTest1
  ${CXX_TEST_PATH}/itkNodeScalarGradientCalculatorTest1
  ${TEMP}/fixedMesh01.vtk
  ${TEMP}/movingMesh01.vtk
  )

ADD_TEST(NodeScalarGradientCalculatorTest2
  ${CXX_TEST_PATH}/itkNodeScalarGradientCalculatorTest2
  )

ADD_TEST(NodeVectorJacobianCalculatorTest1
  ${CXX_TEST_PATH}/itkNodeVectorJacobianCalculatorTest1
  ${TEMP}/fixedMesh01.vtk
  ${TEMP}/movingMesh01.vtk
  )

ADD_TEST(NodeVectorJacobianCalculatorTest2
  ${CXX_TEST_PATH}/itkNodeVectorJacobianCalculatorTest2
  )

ADD_TEST(TriangleBasisSystemTest1
  ${CXX_TEST_PATH}/itkTriangleBasisSystemTest1
 )

ADD_TEST(TriangleListBasisSystemCalculatorTest1
  ${CXX_TEST_PATH}/itkTriangleListBasisSystemCalculatorTest1
 )

ADD_TEST(TriangleBasisSystemCalculatorTest1
  ${CXX_TEST_PATH}/itkTriangleBasisSystemCalculatorTest1
 )

#ADD_TEST(PointLocator2Test1
#  ${CXX_TEST_PATH}/itkPointLocator2Test1
#  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest2.vtk
#  )
 
ADD_TEST(AssignScalarValuesFilterTest1
  ${CXX_TEST_PATH}/itkAssignScalarValuesQuadEdgeMeshFilterTest1
  ${TEMP}/fixedMesh01.vtk 
  ${TEMP}/movingMesh01.vtk 
  ${TEMP}/newMesh01.vtk 
  )

ADD_TEST(RescaleScalarsQuadEdgeMeshFilterTest  
  ${CXX_TEST_PATH}/itkRescaleScalarsQuadEdgeMeshFilterTest
  ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/rescaledmeshWithScalarsGaussianIC4.vtk
  0.5 # min
  1.0 # max
  )

ADD_TEST(HistogramMatchingQuadEdgeMeshFilterTest
  ${CXX_TEST_PATH}/itkHistogramMatchingQuadEdgeMeshFilterTest
  ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk
  ${TEMP}/rescaledmeshWithScalarsGaussianIC4.vtk
  ${TEMP}/histogramMatchedGaussianIC4.vtk
  1024 # NumberOfHistogramLevels
  7 # NumberOfMatchPoints
  )
SET_TESTS_PROPERTIES( HistogramMatchingQuadEdgeMeshFilterTest PROPERTIES DEPENDS RescaleScalarsQuadEdgeMeshFilterTest)
  
ADD_TEST(WarpQuadEdgeMeshFilterTest
  ${CXX_TEST_PATH}/itkWarpQuadEdgeMeshFilterTest
  ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk 
  ${TEMP}/movingMeshWithScalarsGaussianIC4.vtk 
  ${TEMP}/fixedMeshWithDeformationFieldIC4.vtk
  ${TEMP}/warpedFixedMeshWithScalarsGaussianIC4.vtk 
  )
SET_TESTS_PROPERTIES( WarpQuadEdgeMeshFilterTest PROPERTIES DEPENDS ComputeDeformationFieldMeshTest1)
  
ADD_TEST(QuadEdgeMeshSimilarityCalculatorTest
  ${CXX_TEST_PATH}/itkQuadEdgeMeshSimilarityCalculatorTest
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC1.vtk 
  ${TEST_DATA_ROOT}/Input/movingMeshWithScalarsIC1.vtk 
  0 # Label
  0.946329 # Expected Dice for this test
  0.898126 # Expected Jaccard for this test
  )  
ADD_TEST(ScalarQuadEdgeMeshToListAdaptorTest1
  ${CXX_TEST_PATH}/itkScalarQuadEdgeMeshToListAdaptorTest1
  ${TEST_DATA_ROOT}/Input/fixedMeshWithScalarsIC1.vtk 
  10   # NumberOfBins
  0.0  # Min
  10.0 # Max
  )  

IF( FREESURFER_DIR )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest1
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest1
  ${FREESURFER_DIR}/subjects/bert/surf/lh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/lh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest1.vtk
  )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest2
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest2
  ${FREESURFER_DIR}/subjects/bert/surf/lh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/lh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest2.vtk
  )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest3
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest1
  ${FREESURFER_DIR}/subjects/bert/surf/rh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/rh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest3.vtk
  )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest4
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest2
  ${FREESURFER_DIR}/subjects/bert/surf/rh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/rh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest4.vtk
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest2
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest2
  ${TEMP}
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest1
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest1
  ${FREESURFER_DIR}/subjects/bert/surf/lh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/lh.curv.pial
  )
 
  
ENDIF( FREESURFER_DIR )

ADD_TEST(itkSimpleSurface1
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest2
  ${TEST_DATA_ROOT}/Input/SimpleSurface1/lh.sphere
  ${TEST_DATA_ROOT}/Input/SimpleSurface1/lh.sulc
  ${TEMP}/itkSimpleSurface1.vtk
  )

ADD_TEST(itkSimpleSurface2
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest2
  ${TEST_DATA_ROOT}/Input/SimpleSurface2/lh.sphere
  ${TEST_DATA_ROOT}/Input/SimpleSurface2/lh.sulc
  ${TEMP}/itkSimpleSurface2.vtk
  )

ADD_TEST(itkNormalizeMeshTest1
  ${CXX_TEST_PATH}/itkNormalizeScalarsQuadEdgeMeshFilterTest1
  ${TEST_DATA_ROOT}/Input/SimpleSurface1.vtk
  ${TEMP}/SimpleSurfaceNormalized1.vtk
  )


IF(USE_VTK)
  IF(GENERATE_PAPER_SCREENSHOTS)
  
    ADD_TEST(RenderTest1
      ${CXX_TEST_PATH}/MeshRender
      ${TEMP}/meshWithoutScalarsGaussianIC4.vtk
      ${TEMP}/meshWithoutScalarsGaussianIC4.png
      )  
       
    ADD_TEST(ScalarRenderTest1
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/fixedMeshWithScalarsGaussianIC4.vtk
      ${TEMP}/fixedMeshWithScalarsGaussianIC4.png
      )

    ADD_TEST(ScalarRenderTest2
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/movingMeshWithScalarsGaussianIC4.vtk
      ${TEMP}/movingMeshWithScalarsGaussianIC4.png
      )

    ADD_TEST(ScalarRenderTest3
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/registeredMeshWithScalarsGaussianIC4.vtk
      ${TEMP}/registeredMeshWithScalarsGaussianIC4.png
      )
      
    ADD_TEST(ScalarRenderTest4
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/fixedMesh01.vtk
      ${TEMP}/fixedMesh01.png
      )
      
    ADD_TEST(ScalarRenderTest5
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/movingMesh01.vtk
      ${TEMP}/movingMesh01.png
      )
      
    ADD_TEST(ScalarRenderTest6
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/newMesh01.vtk
      ${TEMP}/newMesh01.png
      )
      
    ADD_TEST(ScalarRenderTest7
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/warpedFixedMeshWithScalarsGaussianIC4.vtk
      ${TEMP}/warpedFixedMeshWithScalarsGaussianIC4.png
      )   
   
    ADD_TEST(ScalarRenderTest8
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/rescaledmeshWithScalarsGaussianIC4.vtk
      ${TEMP}/rescaledmeshWithScalarsGaussianIC4.png
      ) 
          
    ADD_TEST(ScalarRenderTest9
      ${CXX_TEST_PATH}/ScalarMeshRender
      ${TEMP}/histogramMatchedGaussianIC4.vtk
      ${TEMP}/histogramMatchedGaussianIC4.png
      ) 
                 
    ADD_TEST(VectorRenderTest1
      ${CXX_TEST_PATH}/VectorMeshRender
      ${TEMP}/meshWithVectorsIC4.vtk
      1   # vector scale factor
      ${TEMP}/meshWithVectorsIC4.png
      )
      
    ADD_TEST(VectorRenderTest2
      ${CXX_TEST_PATH}/VectorMeshRender
      ${TEMP}/fixedMeshWithDeformationFieldIC4.vtk
      1   # vector scale factor
      ${TEMP}/fixedMeshWithDeformationFieldIC4.png
      )
      
  ENDIF(GENERATE_PAPER_SCREENSHOTS)
ENDIF(USE_VTK)

