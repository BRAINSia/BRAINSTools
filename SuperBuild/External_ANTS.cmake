# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

set(proj ANTS)
set(${proj}_GIT_REPOSITORY "git://github.com/BRAINSia/ANTS.git")
set(${proj}_GIT_TAG "55e34f08ef8b2646bd2043498d12ac2773b2b12d")

set(${proj}_DEPENDENCIES ${ITK_EXTERNAL_NAME} SlicerExecutionModel)

ExternalProject_Add(${proj}
  GIT_REPOSITORY ${${proj}_GIT_REPOSITORY}
  GIT_TAG ${${proj}_GIT_TAG}
  SOURCE_DIR ${proj}
  BINARY_DIR ${proj}-build
  "${cmakeversion_external_update}"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
   -Wno-dev
   --no-warn-unused-cli
  ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
  -DUSE_SYSTEM_ITK:BOOL=ON
  -DUSE_SYSTEM_SLICER_EXECUTION_MODEL:BOOL=ON
  -DITK_DIR:PATH=${ITK_DIR}
  ${COMMON_EXTERNAL_PROJECT_ARGS}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -D$ANTS_SUPERBUILD:BOOL=OFF
  ${${proj}_CMAKE_OPTIONS}
  INSTALL_COMMAND ""
  DEPENDS
  ${${proj}_DEPENDENCIES}
  )
ExternalProject_Add_Step(${proj} forcebuild
    COMMAND ${CMAKE_COMMAND} -E remove
    ${CMAKE_CURRENT_BUILD_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1
  )
