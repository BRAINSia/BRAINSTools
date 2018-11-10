#project (GenusZeroImageFilter)

#-- ---------------------------------------------------------------------------
#Output directories.
#
#In ITK4 quad edge mesh classes moved to Algorithms

set(CLP GenusZeroImageFilter)
set(${CLP}_SOURCE vtkImageGenus0MarchingCubes.cxx genus0.cxx)
StandardBRAINSBuildMacro(
  NAME ${CLP} ADDITIONAL_SRCS ${${CLP}_SOURCE} TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES} ${BRAINSSurfaceTools_ITK_LIBRARIES} vtkITK)

set(CLP QuadMeshDecimation)
set(${CLP}_SOURCE QuadMeshDecimation.cxx)
StandardBRAINSBuildMacro(
  NAME ${CLP} TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES} ${BRAINSSurfaceTools_ITK_LIBRARIES} vtkITK)

set(CLP QuadMeshSmoothing)
StandardBRAINSBuildMacro(
  NAME ${CLP} TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES} ${BRAINSSurfaceTools_ITK_LIBRARIES} vtkITK)

set(CLP prepSurface)
StandardBRAINSBuildMacro(
  NAME ${CLP} TARGET_LIBRARIES BRAINSCommonLib ${VTK_LIBRARIES} ${BRAINSSurfaceTools_ITK_LIBRARIES} vtkITK)

if (BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  enable_testing()
  include(CTest)

  set(TESTBIN ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

  set(G0TestName TEST_GenusZeroImageFilter)
  ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME ${G0TestName} COMMAND ${TESTBIN}/GenusZeroImageFilter --inputVolume DATA{${TestData_DIR}/0036195-T1T2_leftTissueClass.nii.gz} --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/${G0TestName}.nii.gz --connectivity 6 --biggestComponent --connectedComponent )

  set(G0CompareTestName COMPARE_GenusZeroImageFilter)
  ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME ${G0CompareTestName} COMMAND ${TESTBIN}/ImageCompare ${CMAKE_CURRENT_BINARY_DIR}/${G0TestName}.nii.gz DATA{${TestData_DIR}/0036195-T1T2_leftTissueClass_genus0.nii.gz} )

  set_tests_properties( ${G0CompareTestName} PROPERTIES DEPENDS ${G0TestName} )

endif ()
