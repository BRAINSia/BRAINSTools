# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

## External_${extProjName}.cmake files can be recurisvely included,
## and cmake variables are global, so when including sub projects it
## is important make the extProjName and proj variables
## appear to stay constant in one of these files.
## Store global variables before overwriting (then restore at end of this file.)
ProjectDependancyPush(CACHED_extProjName ${extProjName})
ProjectDependancyPush(CACHED_proj ${proj})

# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName BRAINSTools) #The find_package known name
set(proj        BRAINSTools) #This local name
set(${extProjName}_REQUIRED_VERSION "")  #If a required version is necessary, then set this, else leave blank

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

# Set dependency list
set(${extProjName}_DEPENDENCIES ITKv4 SlicerExecutionModel ANTs VTK DCMTK JPEG TIFF Boost teem)
#if(${PROJECT_NAME}_BUILD_DICOM_SUPPORT)
#  list(APPEND ${proj}_DEPENDENCIES DCMTK)
#endif()

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT DEFINED ${extProjName}_SOURCE_DIR)
  #message(STATUS "${__indent}Adding project ${proj}")

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  ### --- Project specific additions here
  # message("VTK_DIR: ${VTK_DIR}")
  # message("ITK_DIR: ${ITK_DIR}")
  # message("SlicerExecutionModel_DIR: ${SlicerExecutionModel_DIR}")
  # message("BOOST_INCLUDE_DIR:PATH=${BOOST_INCLUDE_DIR}")

  set(${proj}_CMAKE_OPTIONS
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DUSE_SYSTEM_ITK:BOOL=ON
      -DUSE_SYSTEM_VTK:BOOL=ON
      -DUSE_SYSTEM_DCMTK:BOOL=ON
      -DUSE_SYSTEM_TEEM:BOOL=ON
      -DUSE_SYSTEM_TIFF:BOOL=ON
      -DUSE_SYSTEM_ANTs:BOOL=ON
      -DUSE_SYSTEM_JPEG:BOOL=ON
      -DUSE_SYSTEM_SimpleITK:BOOL=ON
      -DANTs_SOURCE_DIR:PATH=${ANTs_SOURCE_DIR}
      -DANTs_LIBRARY_DIR:PATH=${ANTs_LIBRARY_DIR}
      -DUSE_SYSTEM_BOOST:BOOL=ON
      -DUSE_SYSTEM_SlicerExecutionModel:BOOL=ON
      -DBOOST_DIR:PATH=${BOOST_ROOT}
      -DBOOST_ROOT:PATH=${BOOST_ROOT}
      -DDCMTK_DIR:PATH=${DCMTK_DIR}
      -DDCMTK_config_INCLUDE_DIR:PATH=${DCMTK_DIR}/include
      -DBOOST_ROOT:PATH=${BOOST_ROOT}
      -DBOOST_INCLUDE_DIR:PATH=${BOOST_INCLUDE_DIR}
      -DJPEG_DIR:PATH=${JPEG_DIR}
      -DJPEG_INCLUDE_DIR:PATH=${JPEG_INCLUDE_DIR}
      -DJPEG_LIB_DIR:PATH=${JPEG_LIB_DIR}
      -DJPEG_LIBRARY:PATH=${JPEG_LIBRARY}
      -DJPEG_DIR:PATH=${JPEG_DIR}
      -DJPEG_INCLUDE_DIR:PATH=${JPEG_INCLUDE_DIR}
      -DJPEG_LIBRARY:PATH=${JPEG_LIBRARY}
      -DTIFF_DIR:PATH=${TIFF_DIR}
      -DTIFF_LIBRARY:FILEPATH=${TIFF_LIBRARY}
      -DTIFF_INCLUDE_DIR:PATH=${TIFF_INCLUDE_DIR}
      -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
      -DSuperBuild_BRAINSTools_USE_GIT_PROTOCOL=${${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL}
      -DITK_DIR:PATH=${ITK_DIR}
      -DVTK_DIR:PATH=${VTK_DIR}
      -DTeem_DIR:PATH=${Teem_DIR}
      -DUSE_ANTS:BOOL=ON
      -D${proj}_USE_QT:BOOL=${LOCAL_PROJECT_NAME}_USE_QT
      -DUSE_AutoWorkup:BOOL=ON
      -DUSE_BRAINSABC:BOOL=ON
      -DUSE_BRAINSConstellationDetector:BOOL=ON
      -DUSE_BRAINSContinuousClass:BOOL=ON
      -DUSE_BRAINSCut:BOOL=ON
      -DUSE_BRAINSDemonWarp:BOOL=ON
      -DUSE_BRAINSFit:BOOL=ON
      -DUSE_BRAINSImageConvert:BOOL=ON
      -DUSE_BRAINSInitializedControlPoints:BOOL=ON
      -DUSE_BRAINSLandmarkInitializer:BOOL=ON
      -DUSE_BRAINSMultiModeSegment:BOOL=ON
      -DUSE_BRAINSMush:BOOL=ON
      -DUSE_BRAINSROIAuto:BOOL=ON
      -DUSE_BRAINSResample:BOOL=ON
      -DUSE_BRAINSSnapShotWriter:BOOL=ON
      -DUSE_BRAINSSurfaceTools:BOOL=ON
      -DUSE_BRAINSTransformConvert:BOOL=ON
      -DUSE_DebugImageViewer:BOOL=OFF
      -DUSE_GTRACT:BOOL=ON
      -DUSE_ICCDEF:BOOL=OFF
      -DUSE_ImageCalculator:BOOL=ON
    )

  ### --- End Project specific additions
  set(${proj}_REPOSITORY "${git_protocol}://github.com/BRAINSia/BRAINSTools.git")
  set(${proj}_GIT_TAG "7cfb5b30e35c1419dfc61970e93dd572788f191c")
  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild/ExternalSources/${proj}
    BINARY_DIR ${proj}-build
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -Wno-dev
      --no-warn-unused-cli
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      ${${proj}_CMAKE_OPTIONS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${extProjName}_DEPENDENCIES}
    )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(${extProjName}_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(BRAINSCommonLib_DIR    ${CMAKE_BINARY_DIR}/${proj}-build/BRAINSTools-build/BRAINSCommonLib)
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} ${${extProjName}_REQUIRED_VERSION} REQUIRED)
    if(NOT ${extProjName}_DIR)
      message(FATAL_ERROR "To use the system ${extProjName}, set ${extProjName}_DIR")
    endif()
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
