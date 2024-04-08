# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# ExternalProject_Include_Dependencies
set(extProjName Cppcheck ) #The find_package known name
set(proj        Cppcheck ) #This local name
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
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  ### --- Project specific additions here
  set(CPPCHECK_BUILD_ENVIRONMENT HAVE_RULES=no CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    DESTDIR=${CMAKE_BINARY_DIR}/ PREFIX=Utils )
  ### --- End Project specific additions
  set(${proj}_REPOSITORY ${git_protocol}://github.com/danmar/cppcheck.git)
  set(${proj}_GIT_TAG origin/master)
  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BUILD_IN_SOURCE 1
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CPPCHECK_BUILD_ENVIRONMENT} ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND ${CPPCHECK_BUILD_ENVIRONMENT} ${CMAKE_MAKE_PROGRAM} install
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )
  set(${extProjName}_EXE ${CMAKE_BINARY_DIR}/Utils/bin/cppcheck)
else()
  if(${USE_SYSTEM_${extProjName}})
    find_program(${proj}_EXE cppcheck DOC "Path of Cppcheck program")
    message("USING the system ${extProjName}, set ${extProjName}_EXE=${${extProjName}_EXE}")
  endif()
  # The project is provided using ${extProjName}_EXE, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  ExternalProject_Add_empty(${proj} "${${proj}_DEPENDENCIES}")
endif()

mark_as_superbuild(
  VARS
    ${extProjName}_EXE:PATHPATH
  LABELS
    "FIND_PACKAGE"
)
