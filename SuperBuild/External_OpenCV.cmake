if(DEFINED OpenCV_DIR AND NOT EXISTS ${OpenCV_DIR})
  message(FATAL_ERROR "OpenCV_DIR variable is defined but corresponds to non-existing directory")
endif()


option(${CMAKE_PROJECT_NAME}_USE_NEW_OPENCV "if on, then use the 2011-10-31 version of OpenCV" ON)

if(${CMAKE_PROJECT_NAME}_USE_NEW_OPENCV)
  set(OPENCV_GIT_TAG "FixNeuralNetwork_20111031") # USE THIS FOR UPDATED VERSION
else()
  set(OPENCV_GIT_TAG "339f97dd675af1cee552a6bb4430689c58559c19") # USE THIS FOR 2010-12-14 version
endif()

if(NOT DEFINED OpenCV_DIR)
  set(OpenCV_DEPEND OpenCV)
  set(proj OpenCV)

  ExternalProject_add(${proj}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build

    GIT_REPOSITORY "${git_protocol}://github.com/hjmjohnson/OpenCV.git"
    GIT_TAG ${OPENCV_GIT_TAG}
    CMAKE_ARGS
    --no-warn-unused-cli
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
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
    UPDATE_COMMAND ""
    INSTALL_COMMAND ""
    )
if(${CMAKE_PROJECT_NAME}_USE_NEW_OPENCV)
  set(OpenCV_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
else()
  set(OpenCV_DIR ${CMAKE_BINARY_DIR}/${proj}-build/share/opencv)
endif()
endif(NOT DEFINED OpenCV_DIR)
