
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED SimpleITK_DIR AND NOT EXISTS ${SimpleITK_DIR})
  message(FATAL_ERROR "SimpleITK_DIR variable is defined but corresponds to non-existing directory")
endif()

# Set dependency list
set(SimpleITK_DEPENDENCIES ITKv4 Swig)

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(SimpleITK)

#
#  SimpleITK externalBuild
#
include(ExternalProject)

configure_file(SuperBuild/SimpleITK_install_step.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/SimpleITK_install_step.cmake
  @ONLY)

set(SimpleITK_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/SimpleITK_install_step.cmake)

find_package( PythonInterp REQUIRED )
find_package ( PythonLibs REQUIRED )

#
# On the Helium machine I ran into trouble with
# SimpleITK not being able to find Python.h.
# After sleuthing around I determined that
# PYTHON_INCLUDE_DIRS pointed to the parent of the
# directory containing Python.h, So if that's
# the case I search for it and amend the patch.
if(NOT EXISTS "${PYTHON_INCLUDE_DIRS}/Python.h")
  file(GLOB_RECURSE PYTHON_H "${PYTHON_INCLUDE_DIRS}/Python.h")
  get_filename_component(PYTHON_INCLUDE_DIRS ${PYTHON_H} PATH )
  set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_DIRS})
#  message("PYTHON_INCLUDE_DIRS=${PYTHON_INCLUDE_DIRS}")
endif()

#
# 
ExternalProject_add(SimpleITK
  SOURCE_DIR SimpleITK
  BINARY_DIR SimpleITK-build
  GIT_REPOSITORY ${git_protocol}://itk.org/SimpleITK.git
  GIT_TAG cc7bd86396b3b2da2e54a9465066171d129ae4b4
  "${cmakeversion_external_update}"
  CMAKE_ARGS
    -Wno-dev
#    --no-warn-unused-cli
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
    ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  # SimpleITK does not work with shared libs turned on
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
  -DITK_DIR:PATH=${ITK_DIR}
  -DBUILD_TESTING:BOOL=OFF
  -DBUILD_DOXYGEN:BOOL=OFF
  -DWRAP_PYTHON:BOOL=ON
  -DWRAP_TCL:BOOL=OFF
  -DWRAP_JAVA:BOOL=OFF
  -DWRAP_RUBY:BOOL=OFF
  -DWRAP_LUA:BOOL=OFF
  -DWRAP_CSHARP:BOOL=OFF
  -DWRAP_R:BOOL=OFF
  -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
  -DPYTHON_LIBRARY:STRING=${PYTHON_LIBRARY}
  -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
  -DPYTHON_DEBUG_LIBRARIES:STRING=${PYTHON_DEBUG_LIBRARIES}
  -DSWIG_EXECUTABLE:PATH=${SWIG_EXECUTABLE}
  #
  INSTALL_COMMAND ${SimpleITK_INSTALL_COMMAND}
  #
  DEPENDS ${SimpleITK_DEPENDENCIES}
  )
