if(DEFINED OpenCV_DIR AND NOT EXISTS ${OpenCV_DIR})
  message(FATAL_ERROR "OpenCV_DIR variable is defined but corresponds to non-existing directory")
endif()

if(NOT DEFINED OpenCV_DIR)
  set(OpenCV_DEPEND OpenCV)
  set(proj OpenCV)
  ExternalProject_add(${proj}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build

    GIT_REPOSITORY "${git_protocol}://github.com/hjmjohnson/OpenCV.git"
    GIT_TAG "339f97dd675af1cee552a6bb4430689c58559c19"
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
    UPDATE_COMMAND ""
    )
  set(OpenCV_DIR ${CMAKE_BINARY_DIR}/${proj}-build/share/opencv)
endif(NOT DEFINED OpenCV_DIR)
