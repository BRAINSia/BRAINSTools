add_subdirectory( CompareSurfaces )


set(configure_header_file vtkITKConfigure.h)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${configure_header_file}.in
  ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
  )

include_directories(${CMAKE_CURRENT_BINARY_DIR})

StandardBRAINSBuildMacro(NAME BRAINSApplySurfaceLabels TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

StandardBRAINSBuildMacro(NAME BRAINSMeasureSurface TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

StandardBRAINSBuildMacro(NAME BRAINSSurfaceGeneration
ADDITIONAL_SRCS vtkITKArchetypeImageSeriesReader.cxx vtkITKArchetypeImageSeriesScalarReader.cxx
TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

if (BUILD_TESTING)

set(project_name BRAINSApplySurfaceLabels)

ExternalData_add_test(${project_name}FetchData NAME TEST_BRAINSSurfaceGeneration
  COMMAND $<TARGET_FILE:BRAINSSurfaceGeneration>
          --inputImageFile DATA{${TestData_DIR}/0036195-T1T2_leftTissueClass_genus0.nii.gz}
          --outputSurface 0036195_SurfaceGenerationTest.vtk
          --surfaceValue 0.5
          --genusNumber
          --decimateSurface
          --numberOfElements 5000
          --smoothSurface
          --numIterations 5 )

ExternalData_add_test(${project_name}FetchData NAME COMPARE_BRAINSSurfaceGeneration
  COMMAND $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile 0036195_SurfaceGenerationTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036195_SurfaceGeneration.vtk}
          --vertexLocation
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_BRAINSSurfaceGeneration PROPERTIES DEPENDS TEST_BRAINSSurfaceGeneration )

ExternalData_add_test(${project_name}FetchData NAME TEST_BRAINSApplySurfaceLabels
  COMMAND $<TARGET_FILE:BRAINSApplySurfaceLabels>
          --inputSurface DATA{${TestData_DIR}/0036195_leftSurfaceFeatures.vtk}
          --inputLabelMap DATA{${TestData_DIR}/0036195-T1T2_ACPC_Talairach_LabelMap.nii.gz}
          --outputSurface 0036195_leftTalLabelsTest.vtk
          --cellDataName TalairachLabels )

ExternalData_add_test(${project_name}FetchData  NAME COMPARE_BRAINSApplySurfaceLabels
  COMMAND $<TARGET_FILE:CompareSurfaces>
          --inputSurfaceFile 0036195_leftTalLabelsTest.vtk
          --refSurfaceFile DATA{${TestData_DIR}/0036195_leftTalLabels.vtk}
          --vertexLocation
          --cellData
          --tolerance 1.0e-7 )

set_tests_properties( COMPARE_BRAINSApplySurfaceLabels PROPERTIES DEPENDS TEST_BRAINSApplySurfaceLabels )

ExternalData_add_test(${project_name}FetchData NAME TEST_BRAINSMeasureSurface0
  COMMAND $<TARGET_FILE:BRAINSMeasureSurface>
          --inputSurface DATA{${TestData_DIR}/0036195_leftTalLabels.vtk}
          --writeCsvFile
          --csvFile 0036195_leftTalLabels.csv
          --arrayName TalairachLabels
          --labelNames frontal,temporal,parietal,occipital
          --scanId 0036195
          --subjectId 0036195 )


ExternalData_add_test(${project_name}FetchData NAME TEST_BRAINSMeasureSurface1
  COMMAND $<TARGET_FILE:BRAINSMeasureSurface>
          --inputSurface DATA{${TestData_DIR}/0036195_leftTalLabels.vtk}
          --writeCsvFile
          --csvFile 0036195_leftTalLabels.csv
          --arrayName TalairachLabels
          --labelNames frontal,temporal,parietal,occipital
          --scanId 0036195
          --subjectId 0036195
          --testDepth
          --totalDepthResults 1.96012,2.07585,2.08634,2.55434
          --gyralDepthResults 1.74871,2.0133,1.72714,2.02547
          --sulcalDepthResults 2.28187,2.14289,2.45074,3.3185 )

ExternalData_add_test(${project_name}FetchData  NAME TEST_BRAINSMeasureSurface2
  COMMAND $<TARGET_FILE:BRAINSMeasureSurface>
          --inputSurface DATA{${TestData_DIR}/0036195_leftTalLabels.vtk}
          --writeCsvFile
          --csvFile 0036195_leftTalLabels.csv
          --arrayName TalairachLabels
          --labelNames frontal,temporal,parietal,occipital
          --scanId 0036195
          --subjectId 0036195
          --testArea
          --totalAreaResults 17601.4,10581.9,12902.1,6337.62
          --gyralAreaResults 10622.1,5474.16,6497.48,3745.42
          --sulcalAreaResults 6979.36,5107.7,6404.65,2592.19 )

ExternalData_add_test(${project_name}FetchData NAME TEST_BRAINSMeasureSurface3
  COMMAND $<TARGET_FILE:BRAINSMeasureSurface>
          --inputSurface DATA{${TestData_DIR}/0036195_leftTalLabels.vtk}
          --writeCsvFile
          --csvFile 0036195_leftTalLabels.csv
          --arrayName TalairachLabels
          --labelNames frontal,temporal,parietal,occipital
          --scanId 0036195
          --subjectId 0036195
          --testCurvature
          --totalCurvatureResults 0.0325048,0.0246843,0.0110332,0.0322241
          --gyralCurvatureResults 0.0967193,0.121922,0.100533,0.112673
          --sulcalCurvatureResults -0.0652248,-0.0795294,-0.0797638,-0.0840148 )

ExternalData_Add_Target( ${project_name}FetchData )  # Name of data management target

endif (BUILD_TESTING)

