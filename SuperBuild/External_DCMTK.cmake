
set(proj DCMTK)

# Set dependency list
set(${proj}_DEPENDENCIES "zlib")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(DCMTK_DIR CACHE)
  find_package(DCMTK REQUIRED)
endif()

# Sanity checks
if(DEFINED DCMTK_DIR AND NOT EXISTS ${DCMTK_DIR})
  message(FATAL_ERROR "DCMTK_DIR variable is defined but corresponds to non-existing directory")
endif()

if(NOT DEFINED DCMTK_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  set(EXTERNAL_PROJECT_OPTIONAL_ARGS)

  set(CMAKE_PROJECT_INCLUDE_EXTERNAL_PROJECT_ARG)
  if(CTEST_USE_LAUNCHERS)
    set(CMAKE_PROJECT_INCLUDE_EXTERNAL_PROJECT_ARG
      "-DCMAKE_PROJECT_DCMTK_INCLUDE:FILEPATH=${CMAKE_ROOT}/Modules/CTestUseLaunchers.cmake")
  endif()

  if(UNIX)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DDCMTK_FORCE_FPIC_ON_UNIX:BOOL=ON
      -DDCMTK_WITH_WRAP:BOOL=OFF   # CTK does not build on Mac with this option turned ON due to library dependencies missing
      )
  endif()

  set(${proj}_REPOSITORY ${git_protocol}://github.com/commontk/DCMTK.git)
  set(${proj}_GIT_TAG "f461865d1759854db56e4c840991c81c77e45bb9")

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${DCMTK_REPOSITORY}
    GIT_TAG ${DCMTK_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${proj}-build
    CMAKE_ARGS -Wno-dev --no-warn-unused-cli
    CMAKE_CACHE_ARGS
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/DCMTK-install
      ${CMAKE_PROJECT_INCLUDE_EXTERNAL_PROJECT_ARG}
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DDCMTK_WITH_DOXYGEN:BOOL=OFF
      -DDCMTK_WITH_ZLIB:BOOL=OFF # see CTK github issue #25
      -DDCMTK_WITH_OPENSSL:BOOL=OFF # see CTK github issue #25
      -DDCMTK_WITH_PNG:BOOL=OFF # see CTK github issue #25
      -DDCMTK_WITH_TIFF:BOOL=OFF  # see CTK github issue #25
      -DDCMTK_WITH_XML:BOOL=OFF  # see CTK github issue #25
      -DDCMTK_WITH_ICONV:BOOL=OFF  # see CTK github issue #178
      -DDCMTK_OVERWRITE_WIN32_COMPILER_FLAGS:BOOL=OFF
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )
  set(DCMTK_DIR ${CMAKE_BINARY_DIR}/${proj}-install)

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS DCMTK_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
