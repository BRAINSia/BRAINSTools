project(TextureMeasureFilter)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
       "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
       FORCE)
endif(NOT CMAKE_BUILD_TYPE)

##############################################################
# The Command Line Module

set(CLP TextureMeasureFilter)
set(${CLP}_SOURCE ${CLP}.cxx)
generateclp(${CLP}_SOURCE ${CLP}.xml ${BRAINSLogo})



#$$$$$$$$$  ---  This is going to create a shared library instead of a program
add_library(${CLP}LIB SHARED ${${CLP}_SOURCE})
##PHASE II -#Replace main function with new name called "ModuleEntryPoint" during pre-procssing
set_target_properties (${CLP}LIB PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
target_link_libraries (${CLP}LIB ${ITK_LIBRARIES})

# Build it
add_executable(${CLP} ${${CLP}_SOURCE})
#slicer3_set_plugins_output_path(${CLP})
target_link_libraries(${CLP} ${ITK_LIBRARIES})

#slicer3_install_plugins(${CLP})
#slicer3_install_plugins(${CLPI})

# Test it

#enable_testing()
#include(CTest)
#set(UPDATE_TYPE svn)

#slicer3_add_plugins_test(${CLP}Test1 ${CLP} --help)
#slicer3_add_plugins_test(${CLPI}Test1 ${CLPI} --help)
##############################################################
