project(ImageRegionPlotter)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
       "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
       FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

set(CLP ImageRegionPlotter)
set(${CLP}_SOURCE ${CLP}.cxx)
#generateclp(${CLP}_SOURCE ${CLP}.xml ${BRAINSLogo})

# Build it
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "BRAINSCutCommon")
#add_executable(${CLP} ${${CLP}_SOURCE})
target_link_libraries(${CLP}  BRAINSCommonLib ${ITK_LIBRARIES})

install(TARGETS ImageRegionPlotter
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib/static)

