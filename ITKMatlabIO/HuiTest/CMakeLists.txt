find_package(ITK)
include(${ITK_USE_FILE})

include_directories(${ITK_INCLUDE_DIRS})
add_executable(testReadNrrd testReadNrrd.cpp ConvertImage.cpp)
target_link_libraries(testReadNrrd ${ITK_LIBRARIES} ${VTK_LIBRARIES})
