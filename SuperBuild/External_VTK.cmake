
set(proj VTK)

set(${proj}_REQUIRED_VERSION "7.10")  #If a required version is necessary, then set this, else leave blank
set(VTK_VERSION_MAJOR 7)

# Set dependency list
set(${proj}_DEPENDENCIES "zlib")

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
  set(BUILD_SHARED_LIBS OFF)
  set(VTK_WRAP_PYTHON OFF)

  if(NOT APPLE)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      #-DDESIRED_QT_VERSION:STRING=4 # Unused
      -DVTK_USE_GUISUPPORT:BOOL=ON
      -DVTK_USE_QVTK_QTOPENGL:BOOL=${${PRIMARY_PROJECT_NAME}_USE_QT}
      -DVTK_Group_Qt:BOOL=${${PRIMARY_PROJECT_NAME}_USE_QT} ##VTK6
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      -DVTK_REQUIRED_OBJCXX_FLAGS:STRING="" # Should not be needed, but is always causing problems on mac
                                            # This is to prevent the garbage collection errors from creeping back in
      )
  else()
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DVTK_USE_CARBON:BOOL=OFF
      -DVTK_USE_COCOA:BOOL=ON # Default to Cocoa, VTK/CMakeLists.txt will enable Carbon and disable cocoa if needed
      -DVTK_USE_X:BOOL=OFF
      -DVTK_USE_GUISUPPORT:BOOL=ON
      -DVTK_USE_QVTK_QTOPENGL:BOOL=${${PRIMARY_PROJECT_NAME}_USE_QT}
      -DVTK_Group_Qt:BOOL=${${PRIMARY_PROJECT_NAME}_USE_QT}  ## VTK6
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      )
  endif()

  list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DModule_vtkGUISupportQt:BOOL=${${PRIMARY_PROJECT_NAME}_USE_QT}
      -DModule_vtkGUISupportQtOpenGL:BOOL=${${PRIMARY_PROJECT_NAME}_USE_QT})

  if(VTK_WRAP_TCL)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DTCL_INCLUDE_PATH:PATH=${TCL_INCLUDE_PATH}
      -DTK_INCLUDE_PATH:PATH=${TK_INCLUDE_PATH}
      -DTCL_LIBRARY:FILEPATH=${TCL_LIBRARY}
      -DTK_LIBRARY:FILEPATH=${TK_LIBRARY}
      -DTCL_TCLSH:FILEPATH=${TCL_TCLSH}
      )
  endif()

  set(CUSTOM_BUILD_COMMAND)
  if(CMAKE_GENERATOR MATCHES ".*Makefiles.*")
    # Use $(MAKE) as build command to propagate parallel make option
    set(CUSTOM_BUILD_COMMAND BUILD_COMMAND "$(MAKE)")
    set(make_command_definition -DMAKE_COMMAND=$(MAKE) )
  else()
    set(make_command_definition -DMAKE_COMMAND=${CMAKE_MAKE_PROGRAM})
  endif()

  if(UNIX)
    configure_file(SuperBuild/VTK_build_step.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/VTK_build_step.cmake
      @ONLY)
    set(CUSTOM_BUILD_COMMAND BUILD_COMMAND ${CMAKE_COMMAND}
      ${make_command_definition}
      -P ${CMAKE_CURRENT_BINARY_DIR}/VTK_build_step.cmake)
  endif()

  # Slicer settings
  # set(${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
  #   "github.com/Slicer/VTK.git" CACHE STRING "Repository from which to get VTK" FORCE)
  # set(${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
  #   "c88dfedb277969e5f1f6c5349d8f7898610e75f4" CACHE STRING "VTK git tag to use" FORCE)
  #
  # mark_as_advanced(${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG)
  if(NOT DEFINED git_protocol)
    set(git_protocol "git")
  endif()
  # set(${proj}_GIT_REPOSITORY "${git_protocol}://github.com/Slicer/VTK.git" CACHE STRING "Repository from which to get VTK" FORCE)
  # set(${proj}_GIT_TAG "ea7cdc4e0b399be244e79392c67fed068c33e454")  # VTK 20141221
  set(${proj}_GIT_REPOSITORY "${git_protocol}://vtk.org/VTK.git" CACHE STRING "Repository from which to get VTK" FORCE)
  set(${proj}_GIT_TAG "37930340bfbfe3d6068543c7612859b42cea4008")  # VTK 20160824

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build
    GIT_REPOSITORY "${${proj}_GIT_REPOSITORY}"
    GIT_TAG ${${proj}_GIT_TAG}
    ${CUSTOM_BUILD_COMMAND}
    CMAKE_ARGS -Wno-dev --no-warn-unused-cli
    CMAKE_CACHE_ARGS
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DCMAKE_INCLUDE_DIRECTORIES_BEFORE:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DVTK_USE_PARALLEL:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=${VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=ON
      -DVTK_WRAP_TCL:BOOL=${VTK_WRAP_TCL}
      -DModule_vtkIOXML:BOOL=ON
      -DModule_vtkIOXMLParser:BOOL=ON
      ${VTK_QT_ARGS}
      ${VTK_MAC_ARGS}
      # ZLIB
      -D${proj}_USE_SYSTEM_ZLIB:BOOL=ON
      -DZLIB_ROOT:PATH=${ZLIB_ROOT}
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
      INSTALL_COMMAND ""
    )
  ### --- End Project specific additions
  set(${proj}_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build)

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

mark_as_superbuild(VTK_SOURCE_DIR:PATH)

mark_as_superbuild(
  VARS ${proj}_DIR:PATH VTK_VERSION_MAJOR:STRING
  LABELS "FIND_PACKAGE"
  )
