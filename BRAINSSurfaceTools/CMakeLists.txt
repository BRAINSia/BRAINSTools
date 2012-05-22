# Find VTK.
find_package(VTK REQUIRED)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
endif(VTK_FOUND)

include_directories(${CMAKE_CURRENT_LIST_DIR}/BRAINSSurfaceCommon)

add_subdirectory(BRAINSAssignSurfaceFeatures)
add_subdirectory(BRAINSSurfaceFlattening)
add_subdirectory(BRAINSSurfaceGeneration)
add_subdirectory(BRAINSSurfaceRegister)
#add_subdirectory(BRAINSSurfaceStat)
