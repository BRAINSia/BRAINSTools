FIND_PATH(FREESURFER_DIR
  SetUpFreeSurfer.sh
  PATHS
  $ENV{FREESURFER_HOME}
  DOC "Directory where FreeSurfer is installed"
  )

#######################
#######################
#
#  Define Executables
#
#######################
#######################

ADD_EXECUTABLE(itkTriangleListBasisSystemCalculatorTest1_rigid itkTriangleListBasisSystemCalculatorTest1.cxx )
TARGET_LINK_LIBRARIES(itkTriangleListBasisSystemCalculatorTest1_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkTriangleBasisSystemCalculatorTest1_rigid itkTriangleBasisSystemCalculatorTest1.cxx )
TARGET_LINK_LIBRARIES(itkTriangleBasisSystemCalculatorTest1_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkTriangleBasisSystemTest1_rigid itkTriangleBasisSystemTest1.cxx )
TARGET_LINK_LIBRARIES(itkTriangleBasisSystemTest1_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkRegularSphereQuadEdgeMeshSourceTest1 itkRegularSphereQuadEdgeMeshSourceTest1.cxx )
TARGET_LINK_LIBRARIES(itkRegularSphereQuadEdgeMeshSourceTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkFreeSurferBinarySurfaceReaderTest1_rigid itkFreeSurferBinarySurfaceReaderTest1.cxx )
TARGET_LINK_LIBRARIES(itkFreeSurferBinarySurfaceReaderTest1_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkFreeSurferBinarySurfaceReaderTest2_rigid itkFreeSurferBinarySurfaceReaderTest2.cxx )
TARGET_LINK_LIBRARIES(itkFreeSurferBinarySurfaceReaderTest2_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshFunctionTest1 itkMeshFunctionTest1.cxx )
TARGET_LINK_LIBRARIES(itkMeshFunctionTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkLinearInterpolateMeshFunctionTest1_rigid itkLinearInterpolateMeshFunctionTest1.cxx )
TARGET_LINK_LIBRARIES(itkLinearInterpolateMeshFunctionTest1_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkLinearInterpolateMeshFunctionTest2_rigid itkLinearInterpolateMeshFunctionTest2.cxx )
TARGET_LINK_LIBRARIES(itkLinearInterpolateMeshFunctionTest2_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkLinearInterpolateMeshFunctionTest3 itkLinearInterpolateMeshFunctionTest3.cxx )
TARGET_LINK_LIBRARIES(itkLinearInterpolateMeshFunctionTest3 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkLinearInterpolateMeshFunctionTest4 itkLinearInterpolateMeshFunctionTest4.cxx )
TARGET_LINK_LIBRARIES(itkLinearInterpolateMeshFunctionTest4 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkNearestNeighborInterpolateMeshFunctionTest1 itkNearestNeighborInterpolateMeshFunctionTest1.cxx )
TARGET_LINK_LIBRARIES(itkNearestNeighborInterpolateMeshFunctionTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshKdTreeGeneratorTest itkQuadEdgeMeshKdTreeGeneratorTest.cxx )
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshKdTreeGeneratorTest ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkPointLocator2Test1_rigid itkPointLocator2Test1.cxx )
TARGET_LINK_LIBRARIES(itkPointLocator2Test1_rigid ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshToMeshRegistrationMethodTest1 itkMeshToMeshRegistrationMethodTest1.cxx)
TARGET_LINK_LIBRARIES(itkMeshToMeshRegistrationMethodTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshToMeshRegistrationMethodTest2 itkMeshToMeshRegistrationMethodTest2.cxx)
TARGET_LINK_LIBRARIES(itkMeshToMeshRegistrationMethodTest2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshToMeshRegistrationMethodTest3 itkMeshToMeshRegistrationMethodTest3.cxx)
TARGET_LINK_LIBRARIES(itkMeshToMeshRegistrationMethodTest3 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshToMeshMetricTest1 itkMeshToMeshMetricTest1.cxx)
TARGET_LINK_LIBRARIES(itkMeshToMeshMetricTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeanSquaresMeshToMeshMetricTest1 itkMeanSquaresMeshToMeshMetricTest1.cxx)
TARGET_LINK_LIBRARIES(itkMeanSquaresMeshToMeshMetricTest1 ITKQuadEdgeMesh ITKNumerics ITKIO)

ADD_EXECUTABLE(itkMeanSquaresMeshToMeshMetricTest2 itkMeanSquaresMeshToMeshMetricTest2.cxx)
TARGET_LINK_LIBRARIES(itkMeanSquaresMeshToMeshMetricTest2 ITKQuadEdgeMesh ITKNumerics ITKIO)

ADD_EXECUTABLE(itkMeanSquaresMeshToMeshMetricTest3 itkMeanSquaresMeshToMeshMetricTest3.cxx)
TARGET_LINK_LIBRARIES(itkMeanSquaresMeshToMeshMetricTest3 ITKQuadEdgeMesh ITKNumerics ITKIO)

ADD_EXECUTABLE(itkImageToImageRotationRegistrationTest itkImageToImageRotationRegistrationTest.cxx)
TARGET_LINK_LIBRARIES(itkImageToImageRotationRegistrationTest ITKNumerics ITKIO)

ADD_EXECUTABLE(GenerateRegularSphere_rigid GenerateRegularSphere.cxx)
TARGET_LINK_LIBRARIES(GenerateRegularSphere_rigid ITKQuadEdgeMesh ITKNumerics ITKIO) 


###################
###################
#
#  Define Tests
#
###################
###################

MACRO( GENERATE_SPHERE RADIUS SUBDIVISIONS )

ADD_TEST(GenerateIC${SUBDIVISIONS}
  ${CXX_TEST_PATH}/GenerateRegularSphere_rigid
  ${TEMP}/IC${SUBDIVISIONS}.vtk
  ${RADIUS}         # Sphere radius
  ${SUBDIVISIONS}   # Number of subdivisions
  )

ENDMACRO( GENERATE_SPHERE )


GENERATE_SPHERE( 100.0  0 )
GENERATE_SPHERE( 100.0  1 )
GENERATE_SPHERE( 100.0  2 )
GENERATE_SPHERE( 100.0  3 )
GENERATE_SPHERE( 100.0  4 )
GENERATE_SPHERE( 100.0  5 )
GENERATE_SPHERE( 100.0  6 )
GENERATE_SPHERE( 100.0  7 )
GENERATE_SPHERE( 100.0  8 )
GENERATE_SPHERE( 100.0  9 )



ADD_TEST(itkTriangleListBasisSystemCalculatorTest1_rigid
  ${CXX_TEST_PATH}/itkTriangleListBasisSystemCalculatorTest1_rigid
 )

ADD_TEST(itkTriangleBasisSystemCalculatorTest1_rigid
  ${CXX_TEST_PATH}/itkTriangleBasisSystemCalculatorTest1_rigid
 )

ADD_TEST(itkTriangleBasisSystemTest1_rigid
  ${CXX_TEST_PATH}/itkTriangleBasisSystemTest1_rigid
 )

ADD_TEST(itkRegularSphereQuadEdgeMeshSourceTest1
  ${CXX_TEST_PATH}/itkRegularSphereQuadEdgeMeshSourceTest1
  ${TEMP}/itkRegularSphereQuadEdgeMeshSourceTest1.vtk
 )

ADD_TEST(itkMeanSquaresMeshToMeshMetricTest1
  ${CXX_TEST_PATH}/itkMeanSquaresMeshToMeshMetricTest1
  )

ADD_TEST(itkMeanSquaresMeshToMeshMetricTest2
  ${CXX_TEST_PATH}/itkMeanSquaresMeshToMeshMetricTest2
  )

ADD_TEST(itkMeanSquaresMeshToMeshMetricTest3
  ${CXX_TEST_PATH}/itkMeanSquaresMeshToMeshMetricTest3
  )


IF( FREESURFER_DIR )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest1_rigid
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest1_rigid
  ${FREESURFER_DIR}/subjects/bert/surf/lh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/lh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest1.vtk
  )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest2_rigid
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest2_rigid
  ${FREESURFER_DIR}/subjects/bert/surf/lh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/lh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest2.vtk
  )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest3
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest1_rigid
  ${FREESURFER_DIR}/subjects/bert/surf/rh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/rh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest3.vtk
  )

ADD_TEST(itkFreeSurferBinarySurfaceReaderTest4
  ${CXX_TEST_PATH}/itkFreeSurferBinarySurfaceReaderTest2_rigid
  ${FREESURFER_DIR}/subjects/bert/surf/rh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/rh.curv.pial
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest4.vtk
  )

ADD_TEST(itkMeshFunctionTest1
  ${CXX_TEST_PATH}/itkMeshFunctionTest1
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest1_rigid
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest1_rigid
  ${FREESURFER_DIR}/subjects/bert/surf/lh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/lh.curv.pial
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest2_rigid
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest2_rigid
  ${TEMP}
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest3
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest3
  ${TEMP}
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest4a
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest4
  ${TEMP}/IC4.vtk
  100.0
  ${TEMP}/perturbationFieldIC4.vtk
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest4b
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest4
  ${TEMP}/IC5.vtk
  100.0
  ${TEMP}/perturbationFieldIC5.vtk
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest4c
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest4
  ${TEMP}/IC6.vtk
  100.0
  ${TEMP}/perturbationFieldIC6.vtk
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest4d
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest4
  ${TEMP}/IC7.vtk
  100.0
  ${TEMP}/perturbationFieldIC7.vtk
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest4e
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest4
  ${TEMP}/IC8.vtk
  100.0
  ${TEMP}/perturbationFieldIC8.vtk
  )

ADD_TEST(itkLinearInterpolateMeshFunctionTest4f
  ${CXX_TEST_PATH}/itkLinearInterpolateMeshFunctionTest4
  ${TEMP}/IC9.vtk
  100.0
  ${TEMP}/perturbationFieldIC9.vtk
  )

ADD_TEST(itkNearestNeighborInterpolateMeshFunctionTest1
  ${CXX_TEST_PATH}/itkNearestNeighborInterpolateMeshFunctionTest1
  ${FREESURFER_DIR}/subjects/bert/surf/lh.pial
  ${FREESURFER_DIR}/subjects/bert/surf/lh.curv.pial
  )

ADD_TEST(itkQuadEdgeMeshKdTreeGeneratorTest1
  ${CXX_TEST_PATH}/itkQuadEdgeMeshKdTreeGeneratorTest
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest2.vtk
  ${TEMP}/itkQuadEdgeMeshKdTreeGeneratorTest1.vtk
  )

ADD_TEST(itkPointLocator2Test1_rigid
  ${CXX_TEST_PATH}/itkPointLocator2Test1_rigid
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest2.vtk
  )

ADD_TEST(itkMeshToMeshRegistrationMethodTest1a
  ${CXX_TEST_PATH}/itkMeshToMeshRegistrationMethodTest1
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest1.vtk
  )

ADD_TEST(itkMeshToMeshRegistrationMethodTest2
  ${CXX_TEST_PATH}/itkMeshToMeshRegistrationMethodTest2
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest1.vtk
  )

ADD_TEST(itkMeshToMeshRegistrationMethodTest3
  ${CXX_TEST_PATH}/itkMeshToMeshRegistrationMethodTest3
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest1.vtk
  )

ADD_TEST(itkMeshToMeshMetricTest1a
  ${CXX_TEST_PATH}/itkMeshToMeshMetricTest1
  ${TEMP}/itkFreeSurferBinarySurfaceReaderTest1.vtk
  )


ENDIF( FREESURFER_DIR )

