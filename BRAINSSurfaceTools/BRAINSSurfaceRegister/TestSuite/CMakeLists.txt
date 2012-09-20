add_executable(QuadEdgeMeshSimilarityTest QuadEdgeMeshSimilarityTest.cxx)
target_link_libraries(QuadEdgeMeshSimilarityTest ${ITK_LIBRARIES})

##  Test goes here:
ExternalData_add_test(${PROJECT_NAME}FetchData NAME TEST_SurfaceQuadEdgeMeshSimilarity
  COMMAND $<TARGET_FILE:QuadEdgeMeshSimilarityTest>
  DATA{${TestData_DIR}/0436396_left_ManualLabels.vtk}
  DATA{${TestData_DIR}/0436396_left_BRAINSLabels.vtk}
  0   # LabelValue
  0.903387 # ExpectedDice
  )

## - ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
