project(AllFeatureVectors)
cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)
#include(CPack)

enable_testing()
#include(Dart)
include(CTest)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${AllFeatureVectors_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${AllFeatureVectors_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${AllFeatureVectors_BINARY_DIR}/bin)

#find_package( VTK REQUIRED )
#include(${VTK_USE_FILE})

find_package( ITK REQUIRED )
include(${ITK_USE_FILE})

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

##  In many cases sub-projects depending on SlicerExectuion Model
##  that can be built stand alone are combined in larger packages
##  This logic will include SlicerExectionModel only if it
##  has not already been built.
find_package(GenerateCLP NO_MODULE QUIET)
#if( NOT GenerateCLP_DIR )
#      add_subdirectory(SlicerExecutionModel)
#endif( NOT GenerateCLP_DIR )
find_package(GenerateCLP NO_MODULE REQUIRED)
if(GenerateCLP_DIR)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_DIR)
  message(FATAL_ERROR "Can't build without GenerateCLP. Please set GenerateCLP_DIR")
endif(GenerateCLP_DIR)

add_subdirectory(CombineLabelToMask)
add_subdirectory(CannyEdgeLevelSetFilter)
add_subdirectory(DilateMask)
#add_subdirectory(Opening)
#add_subdirectory(Rescaler)
add_subdirectory(TextureMeasureFilter)
