project( BRAINSContinuousClass)

cmake_minimum_required(VERSION 2.8)

add_subdirectory(liblinear-1.8)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/liblinear-1.8)

set(ALL_PROGS_LIST
    BRAINSContinuousClass
  )
foreach(prog ${ALL_PROGS_LIST})
	StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES liblinear ${ITK_LIBRARIES} )
endforeach()

### Add the Testing Subdirectory.
if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif(BUILD_TESTING)

