set(extProjName "qhull")
if(DEFINED qhull_DIR AND NOT EXISTS ${qhull_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

set(QHULL_GIT_REPO "${git_protocol}://gitorious.org/qhull/qhull.git") # USE THIS FOR UPDATED VERSION
set(QHULL_GIT_TAG "master") # USE THIS FOR UPDATED VERSION

if(NOT DEFINED qhull_DIR)
  set(qhull_DEPEND qhull)
  set(proj qhull)

  ExternalProject_add(${proj}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build

    GIT_REPOSITORY ${QHULL_GIT_REPO}
    GIT_TAG ${QHULL_GIT_TAG}
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_ARGS
      -Wno-dev
      --no-warn-unused-cli
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_TESTS:BOOL=OFF
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/${proj}-install
    )
  set(qhull_DIR ${CMAKE_BINARY_DIR}/${proj}-install/share/qhull/)
endif()
