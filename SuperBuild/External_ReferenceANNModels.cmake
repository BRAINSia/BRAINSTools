#
# ReferenceANNModels is a reference set of prior-probabilities and atlas images.
#

# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED ReferenceANNModels_DIR AND NOT EXISTS ${ReferenceANNModels_DIR})
  message(FATAL_ERROR "ReferenceANNModels_DIR variable is defined but corresponds to non-existing directory")
endif()

if(DEFINED ANNMODELS_NAME AND NOT EXISTS ${ReferenceANNModels_DIR}/${ANNMODELS_NAME})
  message(FATAL_ERROR "ANNMODELS_NAME variable is defined but <ReferenceANNModels_DIR>/<ANNMODELS_NAME> corresponds to non-existing directory")
endif()

# Set dependency list
set(ReferenceANNModels_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(ReferenceANNModels)
set(proj ReferenceANNModels)

if(NOT DEFINED ReferenceANNModels_DIR OR NOT DEFINED ANNMODELS_NAME)
  set(ANNMODELS_VERSION 20111027)
  set(ANNMODELS_URL http://www.psychiatry.uiowa.edu/users/hjohnson/ftp/ANNModels_${ANNMODELS_VERSION}.tar.gz)
  set(ANNMODELS_NAME ANNModels/ANNModels_${ANNMODELS_VERSION})

  ExternalProject_add(${proj}
    URL ${ANNMODELS_URL}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    UPDATE_COMMAND ""
    BUILD_COMMAND ""
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DReferenceANNModels_XML_DIR:PATH=<BINARY_DIR>
      -DANNMODELS_VERSION:STRING=${ANNMODELS_VERSION}
    INSTALL_COMMAND ""
    DEPENDS
      ${ReferenceANNModels_DEPENDENCIES}
    )
  set(ReferenceANNModels_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build)

else()
  # The project is provided using ReferenceANNModels_DIR, nevertheless since other
  # project may depend on ReferenceANNModels, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${ReferenceANNModels_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ReferenceANNModels_DIR:PATH ANNMODELS_NAME:STRING)

