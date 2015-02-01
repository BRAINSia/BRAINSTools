# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
# get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
# if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
#   return()
# endif()
# set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

include(${CMAKE_CURRENT_LIST_DIR}/EP_Autoconf_Utils.cmake)

# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# ExternalProject_Include_Dependencies
set(proj        TIFF) #This local name

#if(${USE_SYSTEM_${proj}})
#  unset(${proj}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR})
  message(FATAL_ERROR "${proj}_DIR variable is defined but corresponds to non-existing directory (${${proj}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES OpenJPEG)

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(NOT ( DEFINED "USE_SYSTEM_${proj}" AND "${USE_SYSTEM_${proj}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT:STRING=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${CMAKE_OSX_DEPLOYMENT_TARGET})
    set(APPLE_CFLAGS " -DHAVE_APPLE_OPENGL_FRAMEWORK")
  endif()

  ### --- Project specific additions here

  AutoConf_FLAGS(${proj}_CFLAGS C "${APPLE_CFLAGS}")
  AutoConf_FLAGS(${proj}_CXXFLAGS CXX "${APPLE_CFLAGS}")

  ### --- End Project specific additions
  set(${proj}_REPOSITORY ${git_protocol}://github.com/BRAINSia/libtiff.git)
  set(${proj}_GIT_TAG 7ddcdf9b098102011e5a7ceafe53711fad93fde7)

  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
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
  set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-install)
  set(${proj}_INCLUDE_DIR ${CMAKE_BINARY_DIR}/${proj}-install/include)
  set(${proj}_LIBRARY ${CMAKE_BINARY_DIR}/${proj}-install/lib/${CMAKE_STATIC_LIBRARY_PREFIX}tiff${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
  if(${USE_SYSTEM_${proj}})
    find_package(${proj} REQUIRED)
    message("USING the system ${proj}, set ${proj}_DIR=${${proj}_DIR}")
  endif()
  # The project is provided using ${proj}_DIR, nevertheless since other
  # project may depend on ${proj}, let's add an 'empty' one
  ExternalProject_Add_Empty(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${proj}_DIR:PATH)
