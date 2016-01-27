# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
# get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
# if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
#   return()
# endif()
# set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

## External_${extProjName}.cmake files can be recurisvely included,
## and cmake variables are global, so when including sub projects it
## is important make the extProjName and proj variables
## appear to stay constant in one of these files.
## Store global variables before overwriting (then restore at end of this file.)
superbuild_stack_push(CACHED_extProjName ${extProjName})
superbuild_stack_push(CACHED_proj ${proj})

# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# ExternalProject_Include_Dependencies
set(extProjName OpenCV) #The find_package known name
set(proj        OpenCV) #This local name
set(${extProjName}_REQUIRED_VERSION "")  #If a required version is necessary, then set this, else leave blank

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  ### --- Project specific additions here
  set(${proj}_CMAKE_OPTIONS
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_NEW_PYTHON_SUPPORT:BOOL=OFF
      -DBUILD_TESTS:BOOL=OFF
      -DWITH_FFMPEG:BOOL=OFF
      -DWITH_JASPER:BOOL=OFF
      -DWITH_OPENEXR:BOOL=OFF
      -DWITH_PVAPI:BOOL=OFF
      -DWITH_JPEG:BOOL=OFF
      -DWITH_TIFF:BOOL=OFF
      -DWITH_PNG:BOOL=OFF
## The following might cause build issues, here for testing
      -DENABLE_SSE:BOOL=ON
      -DENABLE_SSE2:BOOL=ON
      -DENABLE_SSE3:BOOL=ON
      -DENABLE_SSE41:BOOL=ON
      -DENABLE_SSE42:BOOL=ON
      -DENABLE_SSSE3:BOOL=ON
## The follwing tries to get rid of OPENCV build issue
      -DBUILD_opencv_calib3d:BOOL=OFF
      -DBUILD_opencv_contrib:BOOL=OFF
      -DBUILD_opencv_core:BOOL=ON
      -DBUILD_opencv_features2d:BOOL=OFF
      -DBUILD_opencv_flann:BOOL=ON
      -DBUILD_opencv_highgui:BOOL=OFF
      -DBUILD_opencv_imgproc:BOOL=OFF
      -DBUILD_opencv_legacy:BOOL=OFF
      -DBUILD_opencv_ml:BOOL=ON
      -DBUILD_opencv_nonfree:BOOL=OFF
      -DBUILD_opencv_objdetect:BOOL=OFF
      -DBUILD_opencv_photo:BOOL=OFF
      -DBUILD_opencv_python:BOOL=OFF
      -DBUILD_opencv_stitching:BOOL=OFF
      -DBUILD_opencv_ts:BOOL=OFF
      -DBUILD_opencv_video:BOOL=OFF
      -DBUILD_opencv_videostab:BOOL=OFF
      -DBUILD_opencv_world:BOOL=OFF

      -DBUILD_opencv_superres:BOOL=OFF
      -DBUILD_opencv_python2:BOOL=OFF
      -DBUILD_opencv_videoio:BOOL=OFF
      -DBUILD_opencv_java:BOOL=OFF
      -DBUILD_opencv_imgcodec:BOOL=OFF

## Turn off GPU supports
      -DWITH_CUDA:BOOL=OFF
      -DWITH_CUFFT:BOOL=OFF
      -DWITH_OPENCL:BOOL=OFF
      -DWITH_OPENCLAMDFFT:BOOL=OFF
      -DWITH_VTK:BOOL=OFF
      -DBUILD_opencv_matlab:BOOL=OFF
      -DWITH_EIGEN:BOOL=OFF
      -DWITH_GIGEAPI:BOOL=OFF
      -DWITH_GSTREAMER:BOOL=OFF
      -DWITH_GTK:BOOL=OFF
      -DWITH_LIBV4L:BOOL=OFF
      -DWITH_OPENCLAMDBLAS:BOOL=OFF
      -DWITH_V4L:BOOL=OFF
      -DWITH_WEBP:BOOL=OFF

      -DWITH_IPP:BOOL=OFF
      -DBUILD_DOCS:BOOL=OFF
      -DBUILD_FAT_JAVA_LIB:BOOL=OFF
      -DBUILD_PERF_TESTS:BOOL=OFF
      -DBUILD_PACKAGE:BOOL=OFF
      -DBUILD_WITH_DEBUG_INFO:BOOL=OFF



      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/${proj}-install
    )

  ### --- End Project specific additions
  #set(${proj}_REPOSITORY "${git_protocol}://github.com/Itseez/opencv")
  #set(${proj}_GIT_TAG "2.4.9") # USE THIS FOR UPDATED VERSION
  set(${proj}_REPOSITORY "${git_protocol}://github.com/BRAINSia/opencv.git") # USE THIS FOR UPDATED VERSION
  set(${proj}_GIT_TAG 1cd3c6f36467d59028bb8e2d7bd0a92e519db044)  # "20160127
  #set(${proj}_GIT_TAG "20140630_Upstream") # USE THIS FOR UPDATED VERSION for GCC 4.4.7 on RHEL6
  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${proj}-build
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS -Wno-dev --no-warn-unused-cli
    CMAKE_CACHE_ARGS
      ${${proj}_CMAKE_OPTIONS}
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/share/OpenCV/)
else()
  if(${USE_SYSTEM_${extProjName}})
    if(NOT ${${extProjName}_DIR})
      message(FATAL_ERROR "${extProjName}_DIR is required but not set.")
    endif()
    if(NOT IS_DIRECTORY "${${extProjName}_DIR}")
      message(FATAL_ERROR "${extProjName}_DIR is set but doesn't exist: ${${extProjName}_DIR}")
    endif()
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  ExternalProject_Add_Empty(${proj} "${${proj}_DEPENDENCIES}")
endif()

mark_as_superbuild(
  VARS
    ${extProjName}_DIR:FILE
  LABELS
     "FIND_PACKAGE"
)

superbuild_stack_pop(CACHED_extProjName extProjName)
superbuild_stack_pop(CACHED_proj proj)
