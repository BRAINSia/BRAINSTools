#
# ReferenceAtlas is a reference set of prior-probabilities and atlas images.
#

# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Include dependent projects if any
set(extProjName ReferenceAtlas) #The find_package known name
set(proj ${extProjName})              #This local name

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

if(DEFINED ATLAS_NAME AND NOT EXISTS ${${extProjName}_DIR}/${ATLAS_NAME})
  message(FATAL_ERROR "ATLAS_NAME variable is defined but <${extProjName}_DIR>/<ATLAS_NAME> corresponds to non-existing directory (${${extProjName}_DIR}/${ATLAS_NAME})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")

SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT DEFINED ${extProjName}_DIR OR NOT DEFINED ATLAS_NAME)

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  ### --- Project specific additions here
  set(ATLAS_VERSION 20130106)
  # REMOVE set(ATLAS_MIDAS_CODE http://slicer.kitware.com/midas3/download?items=12780)
  # REMOVE set(ATLAS_MD5 ed9e635ef8681f2b0c666aa72e77021c)
  set(ATLAS_SVN_REPOSITORY http://www.nitrc.org/svn/brainstestdata/BRAINSAtlas)
  set(ATLAS_SVN_REVISION 55)


  set(${proj}_CMAKE_OPTIONS
      -DReferenceAtlas_XML_DIR:PATH=<BINARY_DIR>
      -DATLAS_VERSION:STRING=${ATLAS_VERSION}
      )
  ### --- End Project specific additions
  ## Midas version of Atlas:  Atlas_${ATLAS_VERSION}.tar.gz
  # REMOVE set(ATLAS_URL "${ATLAS_MIDAS_CODE}?Atlas_${ATLAS_VERSION}.tar.gz")
                                                                    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                                    #  This is ignored by midas, but allows the filename for download
                                                                    #  to be generated.
                                                                    #  It is a hack that seems to work.
                                                                    #  If the atlas needs to be changed, then the items=#### will 
                                                                    #  need to be determined from the slicer.kitware.com
                                                                    #  web page and filled in appropriately.
  set(ATLAS_NAME Atlas/Atlas_${ATLAS_VERSION})
  ExternalProject_add(${proj}
    # REMOVE URL ${ATLAS_URL}
    # REMOVE URL_MD5 ${ATLAS_MD5}
    SVN_REPOSITORY ${ATLAS_SVN_REPOSITORY}
    SVN_REVISION -r ${ATLAS_SVN_REVISION}
    SVN_USERNAME slicerbot
    SVN_PASSWORD slicer
    SVN_TRUST_CERT
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -Wno-dev
      --no-warn-unused-cli
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      ${${proj}_CMAKE_OPTIONS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    BUILD_COMMAND ""
    )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} REQUIRED)
    if(NOT ${extProjName}_DIR)
      message(FATAL_ERROR "To use the system ${extProjName}, set ${extProjName}_DIR")
    endif()
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}v4, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)
list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ATLAS_NAME:STRING)

