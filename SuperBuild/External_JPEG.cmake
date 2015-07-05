# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# ExternalProject_Include_Dependencies
set(extProjName JPEG) #The find_package known name
set(proj        JPEG) #This local name
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

  AutoConf_FLAGS(${proj}_CFLAGS C "")
  AutoConf_FLAGS(${proj}_CXXFLAGS CXX "")

  ### --- End Project specific additions
  set(${proj}_REPOSITORY ${git_protocol}://github.com/BRAINSia/JPeg9A.git)
  set(${proj}_GIT_TAG 43056c8b105e3b3f5fff0da8c1b48e58a99c1725)  # BRAINSTools_CompilerCleanup
  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${proj}-build
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
    --enable-shared=No
    --enable-static=Yes
    CC=${CMAKE_C_COMPILER}
    CXX=${CMAKE_CXX_COMPILER}
    "CFLAGS=${${proj}_CFLAGS}"
    "CXXFLAGS=${${proj}_CXXFLAGS}"
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    INSTALL_DIR ${proj}-install
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )
  set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-install)
  #set(${proj}_INCLUDE_DIR ${CMAKE_BINARY_DIR}/${proj}-install/include)
  #set(${proj}_LIB_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib)
  #set(${proj}_LIBRARY ${${proj}_LIB_DIR}/libjpeg.a)

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
  VARS
    ${proj}_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
