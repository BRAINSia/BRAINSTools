# Copyright 2015-2016. Hans J. Johnson <hans-johnson@uiowa.edu>
# All rights reserved. Use of this source code is governed by
# a BSD-style license which can be found in the LICENSE file.
# \author Hans J. Johnson <hans-johnson@uiowa.edu>

enable_language(C CXX)

# --- Provide good defaults for searching for packages (i.e. ismrmrd)
set(CMAKE_PREFIX_PATH "")
if(CMAKE_PREFIX_PATH)
    list(APPEND CMAKE_PREFIX_PATH "/usr/local")
endif()
if(EXISTS $ENV{CMAKE_PREFIX_PATH})
    list(APPEND CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
endif()
if(EXISTS $ENV{ISMRMRD_HOME})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ISMRMRD_HOME})
endif()
list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
## -message(STATUS "Looking for packages in : ${CMAKE_PREFIX_PATH}")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

find_package(OpenMP QUIET)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
include(${CMAKE_CURRENT_LIST_DIR}/../../CMake/LinAlgSelector.cmake)

## Compiler flags -- TODO This could be better, see ITK it won't work on windows builds
#--if(${CMAKE_BUILD_TYPE} MATCHES "Release")
#--    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math -O3") # -Wall -Wextra -Wmissing-prototypes")
#--    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -ffast-math -O3") # -Wall -Wextra")
#--endif()

#--HJ find_package(Matlab REQUIRED)
#--HJ include_directories(${Matlab_INCLUDE_DIRS})
#--HJ message(STATUS "XX: ${Matlab_INCLUDE_DIRS}")

##--HJ find_package(Armadillo REQUIRED)
##--HJ # check whether ILP64 MKL should is used
##--HJ add_definitions(-DUSE_ARMADILLO)
##--HJ set(ARMADILLO_BLAS_LONG_LONG FALSE)
##--HJ if(EXISTS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp")
##--HJ   # Read and parse armadillo config.hpp to find out whether BLAS uses long long
##--HJ   file(STRINGS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp" _armadillo_blas_long_long REGEX "// #define ARMA_BLAS_LONG_LONG")
##--HJ   if ( NOT _armadillo_blas_long_long )
##--HJ     set(ARMADILLO_BLAS_LONG_LONG TRUE)
##--HJ     message("Armadillo is found to use long long for BLAS calls")
##--HJ   else ()
##--HJ     message("Armadillo found. Note that ARMADILLO_BLAS_LONG_LONG _must_ be defined in ${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp to link against the MKL ILP64 interface.")
##--HJ   endif ()
##--HJ   unset(_armadillo_blas_long_long)
##--HJ endif()

#-- set(USE_FFTWF TRUE)
#-- set(USE_FFTWD TRUE)
#-- find_package(FFTW REQUIRED)

#-- find_package(ITK REQUIRED)
#-- include(${ITK_USE_FILE})

#-- find_package(VTK REQUIRED)
#-- include(${VTK_USE_FILE})

find_package(RTK REQUIRED)
include(${RTK_USE_FILE})


add_library(SR_support FFTWUpsample.cpp OpWeightedL2.cpp)
target_link_libraries(SR_support ${ITK_LIBRARIES} ${RTK_LIBRARIES} ${LINALG_LIBRARIES} ${VTK_LIBRARIES})

add_executable(TestSR SRmain.cpp)
target_link_libraries(TestSR SR_support)

add_executable(FFTWUpsampleTest FFTWUpsampleTest.cpp)
target_link_libraries(FFTWUpsampleTest SR_support)

if(0)
add_executable(TESTFFTW TestFFTW2D.cpp )
target_link_libraries(TESTFFTW SR_support)

add_executable(GradientTest GradientTest.cpp )
target_link_libraries(GradientTest SR_support)
endif()


if(0)  #-- SimpleITK not yet ready for review.
find_package(SimpleITK REQUIRED)
include(${SimpleITK_USE_FILE})

add_library(SimpleSR_support SimpleOpWeightedL2.cpp)
target_link_libraries(SipmleSR_support SR_support)

add_executable(SimpleSR SimpleSR.cpp)
target_link_libraries(SimpleSR SimpleSR_support ${ITK_LIBRARIES} ${SimpleITK_LIBRARIES})
endif()
