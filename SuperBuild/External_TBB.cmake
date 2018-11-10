set(proj        TBB) #This local name

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Set dependency list
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

set(${proj}_CMAKE_OPTIONS
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
  -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
  -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
  -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
  #-DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  )

set(${proj}_REPOSITORY ${git_protocol}://github.com/01org/tbb.git)
set(${proj}_GIT_TAG 2019_U2)  # 20181110
ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  GIT_REPOSITORY ${${proj}_REPOSITORY}
  GIT_TAG ${${proj}_GIT_TAG}
  SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
  BINARY_DIR tbb_downloaded/${${proj}_GIT_TAG}
  DOWNLOAD_COMMAND  "" #, no download
  CONFIGURE_COMMAND "" #, no config
  BUILD_COMMAND     "" #, no build
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

#TODO:  Will need to wrap configuration files for compilers
#       that have non-default names
#set(TBB_MAKE_ARGS "compiler=${CMAKE_CXX_COMPILER}")

## Following instructions from for source package integration
## https://github.com/01org/tbb/tree/tbb_2018/cmake#source-package-integration
include(${CMAKE_CURRENT_LIST_DIR}/tbb_cmake/TBBGet.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/tbb_cmake/TBBBuild.cmake)
tbb_get(TBB_ROOT TBB_LOCAL_SRC_DIR
        SAVE_TO ${SOURCE_DOWNLOAD_CACHE}/${proj}
        RELEASE_TAG ${${proj}_GIT_TAG}
        SYSTEM_NAME ${CMAKE_SYSTEM_NAME}
        CONFIG_DIR TBB_DIR
        SOURCE_CODE)
tbb_build(TBB_ROOT ${TBB_LOCAL_SRC_DIR}
          CONFIG_DIR TBB_DIR  #Need to set TBB_DIR for the find_package, and to propogate to other packages
          MAKE_ARGS ${TBB_MAKE_ARGS})

find_package(TBB REQUIRED tbb tbbmalloc)

if(NOT EXISTS ${TBB_DIR})
  message(FATAL_ERROR "'TBB_DIR:PATH=${TBB_DIR}' does not exist")
endif()

mark_as_superbuild(
  VARS
    TBB_DIR:PATH
  LABELS
     "FIND_PACKAGE"
)
#message(STATUS "ZZZ:TBB_DIR=${TBB_DIR}:ZZZ")
#message(STATUS "ZZZ:TBB_LOCAL_SRC_DIR=${TBB_LOCAL_SRC_DIR}:ZZZ")
