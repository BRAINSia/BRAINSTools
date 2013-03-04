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
set(extProjName JPEG) #The find_package known name
set(proj        JPEG) #This local name
set(${extProjName}_REQUIRED_VERSION "")  #If a required version is necessary, then set this, else leave blank

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT ( DEFINED "${extProjName}_DIR" OR ( DEFINED "${USE_SYSTEM_${extProjName}}" AND NOT "${USE_SYSTEM_${extProjName}}" ) ) )


  set(${proj}_URL
    "http://www.infai.org/jpeg/files?get=jpegsrc.v09a.tar.gz")
  set(${proj}_MD5 "546c94c044a3f546927c955e31316dd5")

  AutoConf_FLAGS(${proj}_CFLAGS C "")
  AutoConf_FLAGS(${proj}_CXXFLAGS CXX "")

  ExternalProject_Add(${proj}
    URL ${${proj}_URL}
    URL_MD5 ${${proj}_MD5}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    INSTALL_DIR ${proj}-install
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
    --enable-shared=No
    --enable-static=Yes
    CC=${CMAKE_C_COMPILER}
    CXX=${CMAKE_CXX_COMPILER}
    "CFLAGS=${${proj}_CFLAGS}"
    "CXXFLAGS=${${proj}_CXXFLAGS}"
    DEPENDS
    ${${proj}_DEPENDENCIES}
  )
  set(${extProjName}_DIR
    ${CMAKE_BINARY_DIR}/${proj}-install)
  set(${extProjName}_INCLUDE_DIR
    ${CMAKE_BINARY_DIR}/${proj}-install/include)
  set(${extProjName}_LIB_DIR
    ${CMAKE_BINARY_DIR}/${proj}-install/lib)
  set(${extProjName}_LIBRARY
    ${${extProjName}_LIB_DIR}/libjpeg.a)

else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} ${${extProjName}_REQUIRED_VERSION} REQUIRED)
    if(NOT ${extProjName}_DIR)
      message(FATAL_ERROR "To use the system ${extProjName}, set ${extProjName}_DIR")
    endif()
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
