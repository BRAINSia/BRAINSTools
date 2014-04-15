# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName BOOST) #The find_package known name
set(proj        Boost) #This local name
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
#if(${PROJECT_NAME}_BUILD_DICOM_SUPPORT)
#  list(APPEND ${proj}_DEPENDENCIES DCMTK)
#endif()

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
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
  set(Boost_Install_Dir ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install)
  set(Boost_Configure_Script ${CMAKE_CURRENT_LIST_DIR}/External_Boost_configureboost.cmake)
  set(Boost_Build_Script ${CMAKE_CURRENT_LIST_DIR}/External_Boost_buildboost.cmake)

  ### --- End Project specific additions
# SVN is too slow SVN_REPOSITORY http://svn.boost.org/svn/boost/trunk
# SVN is too slow SVN_REVISION -r "82586"
  set(${proj}_URL http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz )
  set(${proj}_MD5 93780777cfbf999a600f62883bd54b17 )
  if(CMAKE_COMPILER_IS_CLANGXX)
    set(CLANG_ARG -DCMAKE_COMPILER_IS_CLANGXX:BOOL=ON)
  endif()
  set(BOOST_SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj})
  ExternalProject_Add(${proj}
    URL ${${proj}_URL}
    URL_MD5 ${${proj}_MD5}
    SOURCE_DIR ${BOOST_SOURCE_DIR}
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
    ${CLANG_ARG}
    -DBUILD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}
    -DBOOST_INSTALL_DIR:PATH=${Boost_Install_Dir}
    -P ${Boost_Configure_Script}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_COMMAND}
    -DBUILD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/Boost
    -DBOOST_INSTALL_DIR:PATH=${Boost_Install_Dir} -P ${Boost_Build_Script}
  )
  set(BOOST_ROOT        ${BOOST_SOURCE_DIR})
  set(BOOST_INCLUDE_DIR ${BOOST_SOURCE_DIR})
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${proj} ${${extProjName}_REQUIRED_VERSION} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  ExternalProject_Add_Empty(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)
