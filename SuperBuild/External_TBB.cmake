set(proj        TBB) #This local name

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Set dependency list
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

# https://github.com/01org/tbb.git
set(${proj}_REPOSITORY ${git_protocol}://github.com/01org/tbb.git)
set(${proj}_GIT_TAG     2019_U9)  # 20191010
ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  GIT_REPOSITORY ${${proj}_REPOSITORY}
  GIT_TAG ${${proj}_GIT_TAG}
  SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
  BINARY_DIR ${SOURCE_DOWNLOAD_CACHE}/${${proj}_GIT_TAG}
  #DOWNLOAD_COMMAND  "" #, no download
  CONFIGURE_COMMAND "" #, no config
  BUILD_COMMAND     ${CMAKE_CURRENT_LIST_DIR}/build_tbb.sh -b ${SOURCE_DOWNLOAD_CACHE}/${proj} -p ${SOURCE_DOWNLOAD_CACHE}/${proj}-install
  INSTALL_COMMAND   "" #, no install
  UPDATE_COMMAND    "" #, no update
  LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
  LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
  LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
  LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
  ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS -Wno-dev --no-warn-unused-cli
  CMAKE_CACHE_ARGS
    ${${proj}_CMAKE_OPTIONS}
  #INSTALL_COMMAND ""
  DEPENDS
  ${${proj}_DEPENDENCIES}
  )

set(TBB_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}-install/lib/cmake/tbb)

#if( APPLE )
#  set( TBB_MIN_VERSION "2019.0") ## Actually 2019.0.11002 is needed for when OSX MIN version < 10.12
#else()
#  set( TBB_MIN_VERSION "2017.0")
#endif()
#find_package(TBB ${TBB_MIN_VERSION} REQUIRED COMPONENTS tbb tbbmalloc NO_MODULE PATHS ${TBB_DIR} )

#if(NOT EXISTS ${TBB_DIR})
#  message(FATAL_ERROR "'TBB_DIR:PATH=${TBB_DIR}' does not exist")
#endif()

mark_as_superbuild(
  VARS
    TBB_DIR:PATH
  LABELS
     "FIND_PACKAGE"
)
#message(STATUS "ZZZ:TBB_DIR=${TBB_DIR}:ZZZ")
#message(STATUS "ZZZ:TBB_LOCAL_SRC_DIR=${TBB_LOCAL_SRC_DIR}:ZZZ")
