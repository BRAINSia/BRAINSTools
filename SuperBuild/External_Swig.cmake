# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

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
set(extProjName Swig) #The find_package known name
set(proj        Swig) #This local name

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

set( TARGET_SWIG_VERSION 2.0.9 )
if(NOT SWIG_DIR)
  if(WIN32)
    # swig.exe available as pre-built binary on Windows:
    ExternalProject_Add(${proj}
      URL http://prdownloads.sourceforge.net/swig/swigwin-${TARGET_SWIG_VERSION}.zip
      URL_MD5 a1dc34766cf599f49e2092f7973c85f4
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${TARGET_SWIG_VERSION}
      ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      )

    set(SWIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${TARGET_SWIG_VERSION}) # ??
    set(SWIG_EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${TARGET_SWIG_VERSION}/swig.exe)
    set(Swig_DEPEND Swig)
  else()
    # Set dependency list
    set(${proj}_DEPENDENCIES "PCRE")

    # Include dependent projects if any
    SlicerMacroCheckExternalProjectDependency(${proj})
    #
    # SWIG
    #

    # swig uses bison find it by cmake and pass it down
    find_package(BISON)
    set(BISON_FLAGS "" CACHE STRING "Flags used by bison")
    mark_as_advanced(BISON_FLAGS)

    # follow the standard EP_PREFIX locations
    set(swig_binary_dir ${CMAKE_CURRENT_BINARY_DIR}/Swig-prefix/src/Swig-build)
    set(swig_source_dir ${CMAKE_CURRENT_BINARY_DIR}/Swig-prefix/src/Swig)
    set(swig_install_dir ${CMAKE_CURRENT_BINARY_DIR}/Swig/install)

    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild/External_Swig_configure_step.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/External_Swig_configure_step.cmake
      @ONLY)
    set ( swig_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/External_Swig_configure_step.cmake )

    ExternalProject_Add(${proj}
      URL http://prdownloads.sourceforge.net/swig/swig-${TARGET_SWIG_VERSION}.tar.gz
      URL_MD5  54d534b14a70badc226129159412ea85
      LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
      LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
      LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
      LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
      ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
      CONFIGURE_COMMAND ${swig_CONFIGURE_COMMAND}
      DEPENDS PCRE
      )

    set(SWIG_DIR ${swig_install_dir}/share/swig/${TARGET_SWIG_VERSION})
    set(SWIG_EXECUTABLE ${swig_install_dir}/bin/swig)
    set(Swig_DEPEND Swig)
  endif()
endif()

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
