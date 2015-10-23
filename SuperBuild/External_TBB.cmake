# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# ExternalProject_Include_Dependencies
set(extProjName TBB) #The find_package known name
set(proj        TBB) #This local name
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

  ### --- End Project specific additions
  set(${proj}_REPOSITORY ${git_protocol}://github.com/intel-tbb/intel-tbb.git)
  set(${proj}_GIT_TAG master)  # BRAINSTools_CompilerCleanup


  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-bld
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/${proj} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} "CFLAGS=${${proj}_CFLAGS}" "CXXFLAGS=${${proj}_CXXFLAGS}" compiler=clang stdlib=libc++ tbb_root=${CMAKE_CURRENT_BINARY_DIR}/${proj} tbb_build_dir=${CMAKE_CURRENT_BINARY_DIR}/${proj}-bld tbb_build_prefix=tbb
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    INSTALL_COMMAND  cp -R ${CMAKE_CURRENT_BINARY_DIR}/${proj}/include/tbb ${CMAKE_CURRENT_BINARY_DIR}/${proj}-bld/tbb_release
    INSTALL_DIR ""
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
