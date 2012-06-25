
add_executable(QuadEdgeMeshSimilarityTest QuadEdgeMeshSimilarityTest.cxx)
target_link_libraries(QuadEdgeMeshSimilarityTest ${ITK_LIBRARIES})

##  Test goes here:
set(QuadEdgeMeshSimilarityTestProg ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/QuadEdgeMeshSimilarityTest)

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceQuadEdgeMeshSimilarity1
  COMMAND ${LAUNCH_EXE} ${QuadEdgeMeshSimilarityTestProg}
          DATA{${TestData_DIR}/0436396_left_ManualLabels.vtk}
          DATA{${TestData_DIR}/0436396_left_BRAINSLabels.vtk}
          0   # LabelValue
          0.948012 # ExpectedDice
          )

ExternalData_add_test(${project_name}FetchData NAME TEST_SurfaceQuadEdgeMeshSimilarity2
  COMMAND ${LAUNCH_EXE} ${QuadEdgeMeshSimilarityTestProg}
          DATA{${TestData_DIR}/0436396_left_ManualLabels.vtk}
          DATA{${TestData_DIR}/0436396_left_BRAINSLabels.vtk}
          1   # LabelValue
          0.822414 # ExpectedDice
          )

MakeTestDriverFromSEMTool(LabelMaps LabelMapsTest.cxx)

## Should provide exactly the same result as LabelMap
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME COMPARE_SurfaceLabelMap
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:LabelMapsTestDriver>
  --compare
  DATA{${TestData_DIR}/LabelMap.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/LabelMapTest.nii.gz
  --compareIntensityTolerance 13
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1
  LabelMapsTest
  --labelFileList DATA{${TestData_DIR}/l_ACiG_man_gm.mask},DATA{${TestData_DIR}/l_AG_man_gm.mask}
  --outputImageFile ${CMAKE_CURRENT_BINARY_DIR}/LabelMapTest.nii.gz
  )

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
