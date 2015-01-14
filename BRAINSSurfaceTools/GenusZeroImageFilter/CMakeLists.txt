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

if( BUILD_TESTING )
  enable_testing()
  include(CTest)

  set(TESTBIN ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

  add_test( TEST_GenusZeroImageFilter ${TESTBIN}/GenusZeroImageFilter --inputVolume /scratch/BRAINSTools/regressiontest/surface/0036195-T1T2_leftTissueClass.nii.gz --outputVolume /scratch/BRAINSTools-build/Testing/results/GenusZerotest.nii.gz --vtkOutput /scratch/BRAINSTools-build/Testing/results/GenusZerotest.vtk --connectivity 6 --biggestComponent --connectedComponent )

  add_test( COMPARE_GenusZeroImageFilter ${TESTBIN}/ImageCompare /scratch/BRAINSTools-build/Testing/results/GenusZerotest.nii.gz /scratch/BRAINSTools/regressiontest/surface/0036195-T1T2_leftTissueClass_genus0.nii.gz )

  set_tests_properties( COMPARE_GenusZeroImageFilter PROPERTIES DEPENDS TEST_GenusZeroImageFilter )

  endif(BUILD_TESTING)
