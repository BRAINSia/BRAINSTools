#
# ReferenceAtlas is a reference set of prior-probabilities and atlas images.
#

# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED ReferenceAtlas_DIR AND NOT EXISTS ${ReferenceAtlas_DIR})
  message(FATAL_ERROR "ReferenceAtlas_DIR variable is defined but corresponds to non-existing directory")
endif()

if(DEFINED ATLAS_NAME AND NOT EXISTS ${ReferenceAtlas_DIR}/${ATLAS_NAME})
  message(FATAL_ERROR "ATLAS_NAME variable is defined but <ReferenceAtlas_DIR>/<ATLAS_NAME> corresponds to non-existing directory")
endif()

# Set dependency list
set(ReferenceAtlas_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(ReferenceAtlas)
set(proj ReferenceAtlas)

if(NOT DEFINED ReferenceAtlas_DIR OR NOT DEFINED ATLAS_NAME)
  set(ATLAS_VERSION 20120104)
  set(ATLAS_URL http://www.psychiatry.uiowa.edu/users/brainstestdata/Atlas_${ATLAS_VERSION}.tar.gz)
  set(ATLAS_NAME Atlas/Atlas_${ATLAS_VERSION})

  ExternalProject_add(${proj}
    URL ${ATLAS_URL}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    UPDATE_COMMAND ""
    BUILD_COMMAND ""
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DReferenceAtlas_XML_DIR:PATH=<BINARY_DIR>
      -DATLAS_VERSION:STRING=${ATLAS_VERSION}
    INSTALL_COMMAND ""
    DEPENDS
      ${ReferenceAtlas_DEPENDENCIES}
    )
  set(ReferenceAtlas_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build)

else()
  # The project is provided using ReferenceAtlas_DIR, nevertheless since other
  # project may depend on ReferenceAtlas, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${ReferenceAtlas_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ReferenceAtlas_DIR:PATH ATLAS_NAME:STRING)

