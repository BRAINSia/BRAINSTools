# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

include(${CMAKE_CURRENT_LIST_DIR}/EP_Autoconf_Utils.cmake)

## External_${extProjName}.cmake files can be recurisvely included,
## and cmake variables are global, so when including sub projects it
## is important make the extProjName and proj variables
## appear to stay constant in one of these files.
## Store global variables before overwriting (then restore at end of this file.)
ProjectDependancyPush(CACHED_extProjName ${extProjName})
ProjectDependancyPush(CACHED_proj ${proj})

# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName TIFF) #The find_package known name
set(proj        TIFF) #This local name

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES OpenJPEG)

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
    set(APPLE_CFLAGS " -DHAVE_APPLE_OPENGL_FRAMEWORK")
  endif()

  ### --- Project specific additions here

  AutoConf_FLAGS(${proj}_CFLAGS C "${APPLE_CFLAGS}")
  AutoConf_FLAGS(${proj}_CXXFLAGS CXX "${APPLE_CFLAGS}")

  ### --- End Project specific additions
  set(${proj}_URL "http://download.osgeo.org/libtiff/tiff-4.0.3.tar.gz")
  set(${proj}_MD5 "051c1068e6a0627f461948c365290410")
  ExternalProject_Add(${proj}
    URL ${${proj}_URL}
    ${URL_HASH_CLAUSE}
    URL_MD5 ${${proj}_MD5}
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/ExternalSources/${proj}
    BINARY_DIR ${proj}-build
    INSTALL_DIR ${proj}-install
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
    --enable-shared=No
    --enable-static=Yes
    --disable-lzma
    --with-jpeg-lib-dir=${JPEG_LIB_DIR}
    --with-jpeg-include-dir=${JPEG_INCLUDE_DIR}
    CC=${CMAKE_C_COMPILER}
    CXX=${CMAKE_CXX_COMPILER}
    "CFLAGS=${${proj}_CFLAGS}"
    "CXXFLAGS=${${proj}_CXXFLAGS}"
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )

  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install)
  set(${extProjName}_INCLUDE_DIR
    ${CMAKE_BINARY_DIR}/${proj}-install/include)
  set(${extProjName}_LIBRARY
    ${CMAKE_BINARY_DIR}/${proj}-install/lib/libtiff.a)
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
