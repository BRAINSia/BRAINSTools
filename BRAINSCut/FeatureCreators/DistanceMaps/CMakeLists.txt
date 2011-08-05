project(DistanceMaps)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE "Debug")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
       "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
       FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# We need Slicer3
#set(Slicer3_DIR "/export/ni/apps/linux-x86_64/Slicer3-build/")
#if(NOT Slicer3_SOURCE_DIR)
#  find_package(Slicer3 REQUIRED)
#  include(${Slicer3_USE_FILE})
#  #slicer3_set_default_install_prefix_for_external_projects()
#endif(NOT Slicer3_SOURCE_DIR)

# Set reasonable default install prefix and output paths
# (after setting Slicer3_DIR, delete CMAKE_INSTALL_PREFIX and re-configure)

#slicer3_set_default_install_prefix_for_external_projects()

##############################################################
# The Command Line Module

set(CLP DistanceMaps)
set(${CLP}_SOURCE ${CLP}.cxx)

generateclp(${CLP}_SOURCE ${CLP}.xml ${BRAINSLogo})

# Build it

add_executable(${CLP} ${${CLP}_SOURCE})
#slicer3_set_plugins_output_path(${CLP})
target_link_libraries(${CLP} ${ITK_LIBRARIES})


# Install it

#slicer3_install_plugins(${CLP})
#slicer3_install_plugins(${CLPI})

# Test it

#enable_testing()
#include(CTest)
#set(UPDATE_TYPE svn)

#slicer3_add_plugins_test(${CLP}Test1 ${CLP} --help)
#slicer3_add_plugins_test(${CLPI}Test1 ${CLPI} --help)
##############################################################
