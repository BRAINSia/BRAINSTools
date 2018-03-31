add_executable(QuadEdgeMeshSimilarityTest QuadEdgeMeshSimilarityTest.cxx)
target_link_libraries(QuadEdgeMeshSimilarityTest ${BRAINSSurfaceTools_ITK_LIBRARIES} ${VTK_LIBRARIES})
set_target_properties(QuadEdgeMeshSimilarityTest PROPERTIES FOLDER ${MODULE_FOLDER})

##  Test goes here:
ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME TEST_SurfaceQuadEdgeMeshSimilarity
  COMMAND $<TARGET_FILE:QuadEdgeMeshSimilarityTest>
  DATA{${TestData_DIR}/0436396_left_ManualLabels.vtk}
  DATA{${TestData_DIR}/0436396_left_BRAINSLabels.vtk}
  0   # LabelValue
  0.903387 # ExpectedDice
  )
