project(AllFeatureVectors)
cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)
#include(CPack)

enable_testing()
include(CTest)

add_subdirectory(CombineLabelToMask)
add_subdirectory(CannyEdgeLevelSetFilter)
add_subdirectory(DilateMask)
