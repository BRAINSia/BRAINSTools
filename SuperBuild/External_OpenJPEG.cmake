
# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName OpenJPEG) #The find_package known name
set(proj        OpenJPEG) #This local name
set(${extProjName}_REQUIRED_VERSION "")  #If a required version is necessary, then set this, else leave blank

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES zlib)

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  ### --- Project specific additions here
  set(${proj}_CMAKE_OPTIONS
      -DTIFF_DIR:PATH=
      -DZLIB_ROOT:PATH=${ZLIB_ROOT}
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
      -DPNG_PNG_INCLUDE_DIR:STRING=
      -DPNG_LIBRARY_DEBUG:STRING=
      -DPNG_LIBRARY_RELEASE:STRING=
    )

  ### --- End Project specific additions
  set(${proj}_REPOSITORY "${git_protocol}://github.com/uclouvain/openjpeg.git")
  set(${proj}_GIT_TAG "0bda7188b7b545232a341f1d978b1e4feda46fc2") # 20210122

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
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS -Wno-dev --no-warn-unused-cli
    CMAKE_CACHE_ARGS
      ${${proj}_CMAKE_OPTIONS}
      ${EXTERNAL_PROJECT_DEFAULTS}
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/openjpeg-2.0)
  set(${extProjName}_INCLUDE_DIR ${CMAKE_BINARY_DIR}/${proj}-install/include/openjpeg-2.0)
  set(${extProjName}_LIB_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib)
  set(${extProjName}_LIBRARY ${${extProjName}_LIB_DIR}/libopenjp2.a)

else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} ${${extProjName}_REQUIRED_VERSION} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  ExternalProject_Add_Empty(${proj} "${${proj}_DEPENDENCIES}")
endif()

mark_as_superbuild(
  VARS
    ${proj}_INCLUDE_DIR:PATH
    ${proj}_LIB_DIR:PATH
    ${proj}_LIBRARY:FILEPATH
  LABELS "FIND_PACKAGE"
  )

ExternalProject_Message(${proj} "${proj}_INCLUDE_DIR:${${proj}_INCLUDE_DIR}")
ExternalProject_Message(${proj} "${proj}_LIB_DIR:${${proj}_LIB_DIR}")
ExternalProject_Message(${proj} "${proj}_LIBRARY:${${proj}_LIBRARY}")
