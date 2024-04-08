
set(proj VTK)

set(VTK_VERSION_MAJOR 9)
set(VTK_VERSION_MINOR 3)
set(${proj}_REQUIRED_VERSION "9.3")  #If a required version is necessary, then set this, else leave blank

# Set dependency list
set(${proj}_DEPENDENCIES "zlib" )
#"TBB")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(VTK_DIR CACHE)
  unset(VTK_SOURCE_DIR CACHE)
  find_package(VTK REQUIRED NO_MODULE)
endif()

# Sanity checks
if(DEFINED VTK_DIR AND NOT EXISTS ${VTK_DIR})
  message(FATAL_ERROR "VTK_DIR variable is defined but corresponds to non-existing directory")
endif()

if(DEFINED VTK_SOURCE_DIR AND NOT EXISTS ${VTK_SOURCE_DIR})
  message(FATAL_ERROR "VTK_SOURCE_DIR variable is defined but corresponds to non-existing directory")
endif()


if((NOT DEFINED VTK_DIR OR NOT DEFINED VTK_SOURCE_DIR) AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(EXTERNAL_PROJECT_OPTIONAL_ARGS)

  set(VTK_WRAP_TCL OFF)
  set(VTK_WRAP_PYTHON OFF)

  if(NOT APPLE)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      #-DDESIRED_QT_VERSION:STRING=4 # Unused
      -DVTK_USE_GUISUPPORT:BOOL=ON
      -DVTK_USE_QVTK_QTOPENGL:BOOL=${${SUPERBUILD_TOPLEVEL_PROJECT}_USE_QT}
      -DVTK_Group_Qt:BOOL=${${SUPERBUILD_TOPLEVEL_PROJECT}_USE_QT} ##VTK6
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      -DVTK_REQUIRED_OBJCXX_FLAGS:STRING= # Should not be needed, but is always causing problems on mac
                                          # This is to prevent the garbage collection errors from creeping back in
      )
  else()
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DVTK_USE_CARBON:BOOL=OFF
      -DVTK_USE_COCOA:BOOL=ON # Default to Cocoa, VTK/CMakeLists.txt will enable Carbon and disable cocoa if needed
      -DVTK_USE_X:BOOL=OFF
      -DVTK_USE_GUISUPPORT:BOOL=ON
      -DVTK_USE_QVTK_QTOPENGL:BOOL=${${SUPERBUILD_TOPLEVEL_PROJECT}_USE_QT}
      -DVTK_Group_Qt:BOOL=${${SUPERBUILD_TOPLEVEL_PROJECT}_USE_QT}  ## VTK6
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      )
  endif()

  list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DModule_vtkGUISupportQt:BOOL=${${SUPERBUILD_TOPLEVEL_PROJECT}_USE_QT}
      -DModule_vtkGUISupportQtOpenGL:BOOL=${${SUPERBUILD_TOPLEVEL_PROJECT}_USE_QT})

  # Slicer settings
  # set(${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
  #   "github.com/Slicer/VTK.git" CACHE STRING "Repository from which to get VTK" FORCE)
  # set(${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
  #   "c88dfedb277969e5f1f6c5349d8f7898610e75f4" CACHE STRING "VTK git tag to use" FORCE)
  #
  # mark_as_advanced(${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG)
  if(NOT DEFINED git_protocol)
    set(git_protocol "https")
  endif()
  set(vtk_git_protocol "https")
  set(${proj}_GIT_REPOSITORY "${vtk_git_protocol}://gitlab.kitware.com/vtk/VTK.git" CACHE STRING "Repository from which to get VTK" FORCE)
  set(${proj}_GIT_TAG "v9.3.0")   # 20231130

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-${EXTERNAL_PROJECT_BUILD_TYPE}-build
    GIT_REPOSITORY "${${proj}_GIT_REPOSITORY}"
    GIT_TAG ${${proj}_GIT_TAG}
    CMAKE_ARGS -Wno-dev --no-warn-unused-cli
    CMAKE_CACHE_ARGS
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
      ${EXTERNAL_PROJECT_DEFAULTS}
      -DCMAKE_CXX_FLAGS:STRING=${ep_CMAKE_CXX_FLAGS}
      -DCMAKE_C_FLAGS:STRING=${ep_CMAKE_C_FLAGS}
      -DCMAKE_INCLUDE_DIRECTORIES_BEFORE:BOOL=OFF
      -DVTK_USE_PARALLEL:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=${VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=ON
      -DVTK_WRAP_TCL:BOOL=OFF
      -DVTK_WRAP_PYTHON:BOOL=${VTK_WRAP_PYTHON}
      -DModule_vtkIOXML:BOOL=ON
      -DModule_vtkIOXMLParser:BOOL=ON
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
      #INSTALL_COMMAND ""
      DEPENDS "${${proj}_DEPENDENCIES}"
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
