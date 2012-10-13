# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

set(proj NIPYPE)
set(${proj}_GIT_REPOSITORY "git://github.com/nipy/nipype.git")
set(${proj}_GIT_TAG "0a36f46861d165540715c7268bd270b861bc0f83") ## 2012-10-13

ExternalProject_Add(${proj}
  GIT_REPOSITORY ${${proj}_GIT_REPOSITORY}
  GIT_TAG ${${proj}_GIT_TAG}
  SOURCE_DIR ${proj}
  "${cmakeversion_external_update}"
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  )
