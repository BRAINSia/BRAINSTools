# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
set(extProjName Uncrustify) #The find_package known name
set(proj        Uncrustify) #This local name
set(${extProjName}_REQUIRED_VERSION "")  #If a required version is necessary, then set this, else leave blank

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_EXE CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_EXE AND NOT EXISTS ${${extProjName}_EXE})
  message(FATAL_ERROR "${extProjName}_EXE variable is defined but corresponds to non-existing file")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  ### --- Project specific additions here
  set(${proj}_CMAKE_OPTIONS
  )

  ### --- End Project specific additions
  set(${proj}_REPOSITORY "${git_protocol}://github.com/bengardner/uncrustify.git")
  set(${proj}_GIT_TAG "60f3681da60462eda539b78e0c6c3eea823481e5")
  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${proj}-${EXTERNAL_PROJECT_BUILD_TYPE}-build
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${proj}-${EXTERNAL_PROJECT_BUILD_TYPE}-build
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${CMAKE_BINARY_DIR}/Utils
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )
  set(${extProjName}_EXE ${CMAKE_BINARY_DIR}/Utils/bin/uncrustify)
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} ${${extProjName}_REQUIRED_VERSION} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_EXE=${${extProjName}_EXE}")
  endif()
  # The project is provided using ${extProjName}_EXE, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

mark_as_superbuild(
  VARS
    ${extProjName}_EXE:FILE
  LABELS
     "FIND_PACKAGE"
)
