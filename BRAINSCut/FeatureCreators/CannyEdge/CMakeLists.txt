project(CannyEdge)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE "Debug")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
       "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
       FORCE)
endif(NOT CMAKE_BUILD_TYPE)

##############################################################
# The Command Line Module

set(CLP CannyEdge)
set(${CLP}_SOURCE ${CLP}.cxx)
generateclp(${CLP}_SOURCE ${CLP}.xml ${BRAINSLogo})

# Build it

add_executable(${CLP} ${${CLP}_SOURCE})
#slicer3_set_plugins_output_path(${CLP})
target_link_libraries(${CLP} ${ITK_LIBRARIES})

