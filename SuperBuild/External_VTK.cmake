
set(proj VTK)

set(VTK_VERSION_MAJOR 9)
set(VTK_VERSION_MINOR 4)
set(${proj}_REQUIRED_VERSION "9.4")  #If a required version is necessary, then set this, else leave blank

# Set dependency list
set(${proj}_DEPENDENCIES "zlib")
if (Slicer_USE_PYTHONQT)
  list(APPEND ${proj}_DEPENDENCIES python)
endif()
if(Slicer_USE_TBB)
  list(APPEND ${proj}_DEPENDENCIES tbb)
endif()

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(VTK_DIR CACHE)
  unset(VTK_SOURCE_DIR CACHE)
  find_package(VTK REQUIRED NO_MODULE)
endif()

# Sanity checks
if(DEFINED VTK_DIR AND NOT EXISTS ${VTK_DIR})
  message(FATAL_ERROR "VTK_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(DEFINED VTK_SOURCE_DIR AND NOT EXISTS ${VTK_SOURCE_DIR})
  message(FATAL_ERROR "VTK_SOURCE_DIR variable is defined but corresponds to nonexistent directory")
endif()


if((NOT DEFINED VTK_DIR OR NOT DEFINED VTK_SOURCE_DIR) AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(EXTERNAL_PROJECT_OPTIONAL_ARGS)

  set(VTK_WRAP_PYTHON OFF)

  if(Slicer_USE_PYTHONQT)
    set(VTK_WRAP_PYTHON ON)
  endif()

  if(Slicer_USE_PYTHONQT)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
      -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
      -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}
      )
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_VTK9_CMAKE_CACHE_ARGS
      -DPython3_ROOT_DIR:PATH=${Python3_ROOT_DIR}
      -DPython3_INCLUDE_DIR:PATH=${Python3_INCLUDE_DIR}
      -DPython3_LIBRARY:FILEPATH=${Python3_LIBRARY}
      -DPython3_LIBRARY_DEBUG:FILEPATH=${Python3_LIBRARY_DEBUG}
      -DPython3_LIBRARY_RELEASE:FILEPATH=${Python3_LIBRARY_RELEASE}
      -DPython3_EXECUTABLE:FILEPATH=${Python3_EXECUTABLE}
      )
  endif()

  # Markups module needs vtkFrenetSerretFrame, which is available in
  # SplineDrivenImageSlicer remote module.
  list(APPEND EXTERNAL_PROJECT_OPTIONAL_VTK9_CMAKE_CACHE_ARGS
    -DVTK_MODULE_ENABLE_VTK_SplineDrivenImageSlicer:BOOL=YES
    )

  list(APPEND EXTERNAL_PROJECT_OPTIONAL_VTK9_CMAKE_CACHE_ARGS
    -DVTK_MODULE_ENABLE_VTK_ChartsCore:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_ViewsContext2D:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingContext2D:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingContextOpenGL2:STRING=YES
    )

  list(APPEND EXTERNAL_PROJECT_OPTIONAL_VTK9_CMAKE_CACHE_ARGS
    -DVTK_MODULE_ENABLE_VTK_GUISupportQt:STRING=YES
    -DVTK_GROUP_ENABLE_Qt:STRING=YES
    )

  list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
    -DVTK_QT_VERSION:STRING=5
    -DVTK_Group_Qt:BOOL=ON
    -DQt5_DIR:FILEPATH=${Qt5_DIR}
    )

  if(Slicer_USE_TBB)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_VTK9_CMAKE_CACHE_ARGS
      -DTBB_DIR:PATH=${TBB_DIR}
      )
  endif()
  if(APPLE)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DVTK_USE_CARBON:BOOL=OFF
      -DVTK_USE_COCOA:BOOL=ON # Default to Cocoa, VTK/CMakeLists.txt will enable Carbon and disable cocoa if needed
      -DVTK_USE_X:BOOL=OFF
      -DVTK_REQUIRED_OBJCXX_FLAGS:STRING=
      #-DVTK_USE_RPATH:BOOL=ON # Unused
      )
  endif()
  if(UNIX AND NOT APPLE)
    find_package(Fontconfig QUIET)
    if(Fontconfig_FOUND)
      list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
        -DModule_vtkRenderingFreeTypeFontConfig:BOOL=ON
        )
    endif()

    # OpenGL_GL_PREFERENCE
    if(NOT "${OpenGL_GL_PREFERENCE}" MATCHES "^(LEGACY|GLVND)$")
      message(FATAL_ERROR "OpenGL_GL_PREFERENCE variable is expected to be set to LEGACY or GLVND")
    endif()
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DOpenGL_GL_PREFERENCE:STRING=${OpenGL_GL_PREFERENCE}
      )
  endif()

  # Disable Tk when Python wrapping is enabled
  if(Slicer_USE_PYTHONQT)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DVTK_USE_TK:BOOL=OFF
      )
  endif()

  list(APPEND EXTERNAL_PROJECT_OPTIONAL_VTK9_CMAKE_CACHE_ARGS
    -DVTK_BUILD_TESTING:STRING=OFF
    -DVTK_MODULE_ENABLE_VTK_AcceleratorsVTKm:BOOL=NO
    -DCMAKE_INSTALL_LIBDIR:STRING=lib # Force value to prevent lib64 from being used on Linux
    -DVTK_MODULE_ENABLE_VTK_GUISupportQtQuick:BOOL=NO
    )

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_REPOSITORY
    "${EP_GIT_PROTOCOL}://github.com/slicer/VTK.git"
    QUIET
    )

  set(_git_tag)
  if("${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}" STREQUAL "9.2")
    set(_git_tag "59ec450206012e86d4855bc669800499254bfc77") # slicer-v9.2.20230607-1ff325c54-2
    set(vtk_dist_info_version "9.2.20230607")
  elseif("${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}" STREQUAL "9.4")
    set(_git_tag "454bb391dff78c6ff463298a5143ab5b4f0aa083") # slicer-v9.4.2-2025-03-26-13acb1a5d
    set(vtk_dist_info_version "9.4.2")
  else()
    message(FATAL_ERROR "error: Unsupported VTK_VERSION_MAJOR.VTK_VERSION_MINOR: ${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}")
  endif()

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_TAG
    ${_git_tag}
    QUIET
    )

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${Slicer_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${Slicer_${proj}_GIT_TAG}"
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-${EXTERNAL_PROJECT_BUILD_TYPE}-build
    CMAKE_ARGS -Wno-dev --no-warn-unused-cli
    CMAKE_CACHE_ARGS
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
      ${EXTERNAL_PROJECT_DEFAULTS}
      -DCMAKE_CXX_FLAGS:STRING=${ep_CMAKE_CXX_FLAGS}
      -DCMAKE_C_FLAGS:STRING=${ep_CMAKE_C_FLAGS}
      -DCMAKE_INCLUDE_DIRECTORIES_BEFORE:BOOL=OFF
      -DVTK_DEBUG_LEAKS:BOOL=${VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=ON
      #-DVTK_USE_RPATH:BOOL=ON # Unused
      -DVTK_WRAP_PYTHON:BOOL=${VTK_WRAP_PYTHON}
      -DModule_vtkIOXML:BOOL=ON
      -DModule_vtkIOXMLParser:BOOL=ON
      -DVTK_USE_PARALLEL:BOOL=ON
      ## Disable vtk dicom that conflicts with itk dcmtk with duplicate symbols
      -DVTK_GROUP_ENABLE_MPI:STRING=DONT_WANT
      -DVTK_MODULE_ENABLE_VTK_DICOMParser:STRING=DONT_WANT
      -DVTK_MODULE_ENABLE_VTK_vtkDICOM:STRING=DONT_WANT
      -DVTK_ENABLE_WRAPPING:BOOL=OFF
      -DVTK_GROUP_ENABLE_Views:STRING=DONT_WANT
      -DVTK_GROUP_ENABLE_Web:STRING=DONT_WANT
      -DVTK_GROUP_ENABLE_Rendering:STRING=DONT_WANT
      ##
      ${VTK_QT_ARGS}
      ${VTK_MAC_ARGS}
      # ZLIB
      -D${proj}_USE_SYSTEM_ZLIB:BOOL=ON
      -DZLIB_ROOT:PATH=${ZLIB_ROOT}
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
      # TBB
      -DVTK_SMP_IMPLEMENTATION_TYPE:STRING=Sequential
      #-DVTK_SMP_IMPLEMENTATION_TYPE:STRING=TBB
      #-DTBB_DIR:PATH=${TBB_DIR}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )
  ### --- End Project specific additions
  set(${proj}_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/vtk-${${proj}_REQUIRED_VERSION})
  #${CMAKE_CURRENT_BINARY_DIR}/${proj}-${EXTERNAL_PROJECT_BUILD_TYPE}-build)

  set(VTK_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})

  set(PNG_INCLUDE_DIR ${VTK_SOURCE_DIR}/Utilities/vtkpng)

  set(PNG_LIBRARY_DIR ${VTK_DIR}/bin)
  if(CMAKE_CONFIGURATION_TYPES)
    set(PNG_LIBRARY_DIR ${PNG_LIBRARY_DIR}/${CMAKE_CFG_INTDIR})
  endif()
  if(WIN32)
    set(PNG_LIBRARY ${PNG_LIBRARY_DIR}/vtkpng.lib)
  elseif(APPLE)
    set(PNG_LIBRARY ${PNG_LIBRARY_DIR}/libvtkpng.dylib)
  else()
    set(PNG_LIBRARY ${PNG_LIBRARY_DIR}/libvtkpng.so)
  endif()

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

#mark_as_superbuild(VARS VTK_SOURCE_DIR:PATH ALL_PROJECTS)

mark_as_superbuild(
  VARS VTK_DIR:PATH
  VTK_VERSION_MAJOR:STRING
  LABELS "FIND_PACKAGE"
  )
