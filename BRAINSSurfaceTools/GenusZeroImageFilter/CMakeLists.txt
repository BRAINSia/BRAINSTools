project(GenusZeroImageFilter)
set(LOCAL_PROJECT_NAME GenusZeroImageFilter)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)
 
enable_testing()
#include(Dart)
include(CTest)
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBuildMacros.cmake)

###
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

###
CHECKIFSLICER3BUILD()  ## Call the convenience macro

if(NOT ITK_FOUND)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

#-----------------------------------------------------------------------------
# Output directories.
#
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBRAINS3BuildMacros.cmake)

add_subdirectory(vtkITK)

include_directories(${GenusZeroImageFilter_SOURCE_DIR}/vtkITK ${GenusZeroImageFilter_BINARY_DIR}/vtkITK)
include_directories(${GenusZeroImageFilter_SOURCE_DIR} ${GenusZeroImageFilter_BINARY_DIR})

# In ITK4 quad edge mesh classes moved to Algorithms
set(MY_ITK_LIBS ${ITK_LIBRARIES})

if(ITK_VERSION_MAJOR LESS "4")
  set(MY_ITK_LIBS ${MY_ITK_LIBS} ITKQuadEdgeMesh)
endif()

set(CLP GenusZeroImageFilter)
set(${CLP}_SOURCE GenusZeroImageFilter.cxx vtkImageGenus0MarchingCubes.cxx genus0.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries(GenusZeroImageFilter
  vtkGraphics vtkCommon vtkIO vtkITK ${MY_ITK_LIBS})

set(CLP QuadMeshDecimation)
set(${CLP}_SOURCE QuadMeshDecimation.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries(QuadMeshDecimation ${MY_ITK_LIBS})

set(CLP QuadMeshSmoothing)
set(${CLP}_SOURCE QuadMeshSmoothing.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries(QuadMeshSmoothing ${MY_ITK_LIBS})

set(CLP prepSurface)
set(${CLP}_SOURCE prepSurface.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries(prepSurface
  iplKernel
  iplObject
  iplPixelTable
  iplUtils
  ${MY_ITK_LIBS}
)

if (BUILD_TESTING)

set(TESTBIN ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_test( TEST_GenusZeroImageFilter ${TESTBIN}/GenusZeroImageFilter --inputVolume
${BrainsRegressionData_DIR}/trunk/surface/0036195-T1T2_leftTissueClass.nii.gz
--outputVolume ${BRAINS_BINARY_DIR}/Testing/results/GenusZerotest.nii.gz --vtkOutput
${BRAINS_BINARY_DIR}/Testing/results/GenusZerotest.vtk --connectivity 6
--biggestComponent --connectedComponent ) 

add_test( COMPARE_GenusZeroImageFilter ${TESTBIN}/ImageCompare
${BRAINS_BINARY_DIR}/Testing/results/GenusZerotest.nii.gz
${BrainsRegressionData_DIR}/trunk/surface/0036195_leftTissueClass_genus0.nii.gz
)

set_tests_properties( COMPARE_GenusZeroImageFilter PROPERTIES DEPENDS TEST_GenusZeroImageFilter )

endif (BUILD_TESTING)
