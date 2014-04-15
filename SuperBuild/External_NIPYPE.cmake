# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# ExternalProject_Include_Dependencies
set(proj NIPYPE) #The find_package known name
set(proj        NIPYPE) #This local name
set(${proj}_REQUIRED_VERSION "")  #If a required version is necessary, then set this, else leave blank

#if(${USE_SYSTEM_${proj}})
#  unset(${proj}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR})
  message(FATAL_ERROR "${proj}_DIR variable is defined but corresponds to non-existing directory (${${proj}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")

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
  endif()

  ### --- Project specific additions here

  set(${proj}_CMAKE_OPTIONS
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
    )

  ### --- End Project specific additions
  set(${proj}_REPOSITORY "${git_protocol}://github.com/nipy/nipype.git")
  set(${proj}_GIT_TAG "cba8e4128336a647cd3fde2fccbed2abc7be4d38")
  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BUILD_COMMAND ""
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )
  set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj})
else()
  if(${USE_SYSTEM_${proj}})
    find_package(${proj} ${${proj}_REQUIRED_VERSION} REQUIRED)
    message("USING the system ${proj}, set ${proj}_DIR=${${proj}_DIR}")
  endif()
  # The project is provided using ${proj}_DIR, nevertheless since other
  # project may depend on ${proj}, let's add an 'empty' one
  ExternalProject_Add_Empty(${proj} "${${proj}_DEPENDENCIES}")
endif()

mark_as_superbuild(
  VARS ${proj}_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
