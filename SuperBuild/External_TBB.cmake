set(proj        TBB) #This local name

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Set dependency list
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(TBB_COMPILERID "gcc")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(TBB_COMPILERID "clang")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(TBB_COMPILERID "icc")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(TBB_COMPILERID "cl")
else()
  #This is probably invalid
  set(TBB_COMPILERID "${CMAKE_CXX_COMPILER_ID}")
endif()


# Outdated by intel https://github.com/01org/tbb.git
# New verions https://github.com/oneapi-src/oneTBB
# HJ private fixes git@github.com:hjmjohnson/oneTBB.git
set(${proj}_REPOSITORY ${git_protocol}://github.com/hjmjohnson/oneTBB.git)
set(${proj}_GIT_TAG     2020_U2_20200411)  # 2020-04-01

ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  GIT_REPOSITORY ${${proj}_REPOSITORY}
  GIT_TAG ${${proj}_GIT_TAG}
  SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
  BINARY_DIR ${SOURCE_DOWNLOAD_CACHE}/${${proj}_GIT_TAG}
  #DOWNLOAD_COMMAND  "" #, no download
  CONFIGURE_COMMAND "" #, no config
  BUILD_COMMAND
    ${CMAKE_CURRENT_LIST_DIR}/build_tbb.sh
       -c "${CMAKE_C_COMPILER}"
       -x "${CMAKE_CXX_COMPILER}"
       -i "${TBB_COMPILERID}"
       -b "${SOURCE_DOWNLOAD_CACHE}/${proj}"
       -p "${CMAKE_INSTALL_PREFIX}"
       -d " -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER} -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS} -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD} -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED} -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS} "
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
    ${EXTERNAL_PROJECT_DEFAULTS}
  DEPENDS
  ${${proj}_DEPENDENCIES}
  )

#set(TBB_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}-install/lib/cmake/tbb)
set(TBB_DIR ${CMAKE_INSTALL_PREFIX}cmake/tbb)

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
