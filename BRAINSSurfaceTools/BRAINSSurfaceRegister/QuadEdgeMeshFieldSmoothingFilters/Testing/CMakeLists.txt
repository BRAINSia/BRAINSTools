SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
SET(TEST_DATA_ROOT ${CMAKE_SOURCE_DIR}/Data)
SET(TEMP ${CMAKE_BINARY_DIR}/Testing/Temporary)

SET(QuadEdgeMeshRigidRegistration ${PROJECT_SOURCE_DIR}/../QuadEdgeMeshRigidRegistration)

FIND_PATH(FREESURFER_DIR
  SetUpFreeSurfer.sh
  PATHS
  $ENV{FREESURFER_HOME}
  DOC "Directory where FreeSurfer is installed"
  )

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/Testing
  )

#######################
#######################
#
#  Define Executables
#
#######################
#######################

ADD_EXECUTABLE(itkMeshGeneratorHelperTest1_smooth itkMeshGeneratorHelperTest1.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest1_smooth ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest2_smooth itkMeshGeneratorHelperTest2.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest2_smooth ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest3_smooth itkMeshGeneratorHelperTest3.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest3_smooth ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkMeshGeneratorHelperTest4_smooth itkMeshGeneratorHelperTest4.cxx )
TARGET_LINK_LIBRARIES(itkMeshGeneratorHelperTest4_smooth ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshVTKPolyDataReader1 itkQuadEdgeMeshVTKPolyDataReader1.cxx )
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshVTKPolyDataReader1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshVTKPolyDataReader2 itkQuadEdgeMeshVTKPolyDataReader2.cxx )
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshVTKPolyDataReader2 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshVectorPixelValuesSmoothingFilterTest1 itkQuadEdgeMeshVectorPixelValuesSmoothingFilterTest1.cxx)
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshVectorPixelValuesSmoothingFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

ADD_EXECUTABLE(itkQuadEdgeMeshScalarPixelValuesSmoothingFilterTest1 itkQuadEdgeMeshScalarPixelValuesSmoothingFilterTest1.cxx)
TARGET_LINK_LIBRARIES(itkQuadEdgeMeshScalarPixelValuesSmoothingFilterTest1 ITKQuadEdgeMesh ITKIO ITKNumerics)

IF(USE_VTK)

  ADD_EXECUTABLE(ScalarMeshRender  ScalarMeshRender.cxx)
  TARGET_LINK_LIBRARIES(ScalarMeshRender 
    vtkRendering vtkFiltering vtkIO vtkHybrid
    ) 

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

ADD_TEST(itkMeshGeneratorHelperTest1_smooth
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest1_smooth
  ${TEMP}/fixedMesh01.vtk
  ${TEMP}/movingMesh01.vtk
  )

ADD_TEST(itkMeshGeneratorHelperTest2_smooth
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest2_smooth
  ${TEMP}/meshWithVectors.vtk
  )

ADD_TEST(itkMeshGeneratorHelperTest3a
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest3_smooth
  ${TEMP}/meshWithVectorsIC1.vtk
  1      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(itkMeshGeneratorHelperTest3b
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest3_smooth
  ${TEMP}/meshWithVectorsIC4.vtk
  4      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(itkMeshGeneratorHelperTest3c
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest3_smooth
  ${TEMP}/meshWithVectorsIC7.vtk
  7      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(itkMeshGeneratorHelperTest4a
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest4_smooth
  ${TEMP}/meshWithScalarsIC1.vtk
  1      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(itkMeshGeneratorHelperTest4b
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest4_smooth
  ${TEMP}/meshWithScalarsIC4.vtk
  4      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(itkMeshGeneratorHelperTest4c
  ${CXX_TEST_PATH}/itkMeshGeneratorHelperTest4_smooth
  ${TEMP}/meshWithScalarsIC7.vtk
  7      # resolution
  100.0  # Scale (radius)
  )

ADD_TEST(itkQuadEdgeMeshVTKPolyDataReader1
  ${CXX_TEST_PATH}/itkQuadEdgeMeshVTKPolyDataReader1
  ${TEMP}/fixedMesh01.vtk
  ${TEMP}/fixedMesh01verification.vtk
  )

ADD_TEST(itkQuadEdgeMeshVTKPolyDataReader2
  ${CXX_TEST_PATH}/itkQuadEdgeMeshVTKPolyDataReader2
  ${TEMP}/meshWithVectors.vtk
  ${TEMP}/meshWithVectorsVerification.vtk
  )

ADD_TEST(ScalarSmoothingFilterTest1
  ${CXX_TEST_PATH}/itkQuadEdgeMeshScalarPixelValuesSmoothingFilterTest1
  ${TEMP}/fixedMesh01.vtk
  ${TEMP}/fixedMesh01Smoothed.vtk
  1.0   # Lambda
  10    # Iterations
  )

ADD_TEST(VectorSmoothingFilterTest1
  ${CXX_TEST_PATH}/itkQuadEdgeMeshVectorPixelValuesSmoothingFilterTest1
  ${TEMP}/meshWithVectors.vtk
  ${TEMP}/meshWithVectorsSmoothed.vtk
  1.0   # Lambda
  10    # Iterations
  )

ADD_TEST(VectorSmoothingFilterTest2
  ${CXX_TEST_PATH}/itkQuadEdgeMeshVectorPixelValuesSmoothingFilterTest1
  ${TEMP}/meshWithVectorsIC1.vtk
  ${TEMP}/meshWithVectorsIC1Smoothed.vtk
  1.0   # Lambda
  10    # Iterations
  )

ADD_TEST(VectorSmoothingFilterTest3
  ${CXX_TEST_PATH}/itkQuadEdgeMeshVectorPixelValuesSmoothingFilterTest1
  ${TEMP}/meshWithVectorsIC4.vtk
  ${TEMP}/meshWithVectorsIC4Smoothed.vtk
  1.0   # Lambda
  10    # Iterations
  )

ADD_TEST(VectorSmoothingFilterTest4
  ${CXX_TEST_PATH}/itkQuadEdgeMeshVectorPixelValuesSmoothingFilterTest1
  ${TEMP}/meshWithVectorsIC7.vtk
  ${TEMP}/meshWithVectorsIC7Smoothed.vtk
  1.0   # Lambda
  10    # Iterations
  )


MACRO(SCALAR_SMOOTH_AND_RENDER INPUTID OUTPUTSUFIX LAMBDA NUMITER)

  ADD_TEST(ScalarSmooth_${INPUTID}${OUTPUTSUFIX}
    ${CXX_TEST_PATH}/SmoothingMeshWithScalars_smooth
    ${TEMP}/meshWithScalars${INPUTID}.vtk
    ${TEMP}/meshWithScalars${INPUTID}${OUTPUTSUFIX}.vtk
    ${LAMBDA}
    ${NUMITER}
    )

  ADD_TEST(ScalarRender_${INPUTID}${OUTPUTSUFIX}
    ${CXX_TEST_PATH}/ScalarMeshRender
    ${TEMP}/meshWithScalars${INPUTID}${OUTPUTSUFIX}.vtk
    ${TEMP}/meshWithScalars${INPUTID}${OUTPUTSUFIX}.png
    )

ENDMACRO(SCALAR_SMOOTH_AND_RENDER)

MACRO(SCALAR_SMOOTH_AND_RENDER_LAMBDASET INPUTID NUMITER)
    SCALAR_SMOOTH_AND_RENDER( ${INPUTID}  _02_${NUMITER} 0.2  ${NUMITER} )
    SCALAR_SMOOTH_AND_RENDER( ${INPUTID}  _04_${NUMITER} 0.4  ${NUMITER} )
    SCALAR_SMOOTH_AND_RENDER( ${INPUTID}  _08_${NUMITER} 0.8  ${NUMITER} )
    SCALAR_SMOOTH_AND_RENDER( ${INPUTID}  _10_${NUMITER} 1.0  ${NUMITER} )
    SCALAR_SMOOTH_AND_RENDER( ${INPUTID}  _20_${NUMITER} 2.0  ${NUMITER} )
    SCALAR_SMOOTH_AND_RENDER( ${INPUTID}  _40_${NUMITER} 4.0  ${NUMITER} )
ENDMACRO(SCALAR_SMOOTH_AND_RENDER_LAMBDASET)


MACRO(VECTOR_SMOOTH INPUTID OUTPUTSUFIX LAMBDA NUMITER )
  ADD_TEST(VectorSmooth_${INPUTID}${OUTPUTSUFIX}
    ${CXX_TEST_PATH}/SmoothingMeshWithVectors
    ${TEMP}/meshWithVectors${INPUTID}.vtk
    ${TEMP}/meshWithVectors${INPUTID}${OUTPUTSUFIX}.vtk
    ${LAMBDA}
    ${NUMITER}
    )
ENDMACRO(VECTOR_SMOOTH)


MACRO(VECTOR_RENDER INPUTID OUTPUTSUFIX VECTOR_SCALE_FACTOR)
  ADD_TEST(VectorRender_${INPUTID}${OUTPUTSUFIX}
    ${CXX_TEST_PATH}/VectorMeshRender
    ${TEMP}/meshWithVectors${INPUTID}${OUTPUTSUFIX}.vtk
    ${VECTOR_SCALE_FACTOR}
    ${TEMP}/meshWithVectors${INPUTID}${OUTPUTSUFIX}.png
    )
ENDMACRO(VECTOR_RENDER)


MACRO(VECTOR_SMOOTH_AND_RENDER INPUTID OUTPUTSUFIX LAMBDA NUMITER VECTOR_SCALE_FACTOR)
  VECTOR_SMOOTH(${INPUTID} ${OUTPUTSUFIX} ${LAMBDA} ${NUMITER})
  VECTOR_RENDER(${INPUTID} ${OUTPUTSUFIX} ${VECTOR_SCALE_FACTOR})
ENDMACRO(VECTOR_SMOOTH_AND_RENDER)


MACRO(VECTOR_SMOOTH_AND_RENDER_LAMBDASET INPUTID NUMITER)
    VECTOR_SMOOTH_AND_RENDER( ${INPUTID}  _02_${NUMITER} 0.2  ${NUMITER} 20 )
    VECTOR_SMOOTH_AND_RENDER( ${INPUTID}  _04_${NUMITER} 0.4  ${NUMITER} 20 )
    VECTOR_SMOOTH_AND_RENDER( ${INPUTID}  _08_${NUMITER} 0.8  ${NUMITER} 20 )
    VECTOR_SMOOTH_AND_RENDER( ${INPUTID}  _10_${NUMITER} 1.0  ${NUMITER} 20 )
    VECTOR_SMOOTH_AND_RENDER( ${INPUTID}  _20_${NUMITER} 2.0  ${NUMITER} 20 )
    VECTOR_SMOOTH_AND_RENDER( ${INPUTID}  _40_${NUMITER} 4.0  ${NUMITER} 20 )
ENDMACRO(VECTOR_SMOOTH_AND_RENDER_LAMBDASET)



IF(GENERATE_PAPER_SCREENSHOTS)
  IF(USE_VTK)

  SCALAR_SMOOTH_AND_RENDER_LAMBDASET( IC1 1 )
  SCALAR_SMOOTH_AND_RENDER_LAMBDASET( IC1 2 )
  SCALAR_SMOOTH_AND_RENDER_LAMBDASET( IC1 5 )

  SCALAR_SMOOTH_AND_RENDER_LAMBDASET( IC4 1 )
  SCALAR_SMOOTH_AND_RENDER_LAMBDASET( IC4 2 )
  SCALAR_SMOOTH_AND_RENDER_LAMBDASET( IC4 5 )

  VECTOR_SMOOTH_AND_RENDER_LAMBDASET( IC1 1 )
  VECTOR_SMOOTH_AND_RENDER_LAMBDASET( IC1 2 )
  VECTOR_SMOOTH_AND_RENDER_LAMBDASET( IC1 5 )

  VECTOR_SMOOTH_AND_RENDER_LAMBDASET( IC4 1 )
  VECTOR_SMOOTH_AND_RENDER_LAMBDASET( IC4 2 )
  VECTOR_SMOOTH_AND_RENDER_LAMBDASET( IC4 5 )

  VECTOR_SMOOTH_AND_RENDER( IC4  _10_20  1.0 20 400 )
  VECTOR_SMOOTH_AND_RENDER( IC4  _10_40  1.0 40 400 )
  VECTOR_SMOOTH_AND_RENDER( IC4  _10_80  1.0 80 400 )

  ENDIF(USE_VTK)
ENDIF(GENERATE_PAPER_SCREENSHOTS)

