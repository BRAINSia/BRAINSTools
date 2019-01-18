
set(proj ITKv5)

# Set dependency list
set(${proj}_DEPENDENCIES "zlib" )
if( NOT USE_SYSTEM_TBB )
  list(APPEND ${proj}_DEPENDENCIES "TBB" )
endif()
if(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)
  list(APPEND ${proj}_DEPENDENCIES "VTK")
endif()
#if(${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT)
  list(APPEND ${proj}_DEPENDENCIES DCMTK)
#endif()

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(ITK_DIR CACHE)
  find_package(ITK 5.0 COMPONENTS ${${CMAKE_PROJECT_NAME}_ITK_COMPONENTS} REQUIRED NO_MODULE)
endif()

# Sanity checks
if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED ITK_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(ITK_VTK_OPTIONS )

  if(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)
    list(APPEND ITK_VTK_OPTIONS
      -DModule_ITKVtkGlue:BOOL=ON
      -DVTK_DIR:PATH=${VTK_DIR}
      )
  endif()

  if(NOT DEFINED git_protocol)
      set(git_protocol "git")
  endif()

  #set(${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY ${git_protocol}://github.com/InsightSoftwareConsortium/ITK.git
  set(${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY ${git_protocol}://github.com/hjmjohnson/ITK.git)
  set(${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG a1a092c0612e00e94518a851ebd259c72e119474 ) #20180408 - FixSpatialObjectsTesting

  set(EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS)

  if(NOT ${CMAKE_PROJECT_NAME}ITKV3_COMPATIBILITY AND CMAKE_CL_64)
    # enables using long long type for indexes and size on platforms
    # where long is only 32-bits (msvc)
    set(EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DITK_USE_64BITS_IDS:BOOL=ON
      )
  endif()

  if(Slicer_USE_PYTHONQT OR Slicer_BUILD_ITKPython)
    # XXX Ensure python executable used for ITKModuleHeaderTest
    #     is the same as Slicer.
    #     This will keep the sanity check implemented in SlicerConfig.cmake
    #     quiet.
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
      )
  endif()

  if(Slicer_BUILD_ITKPython)

    # Sanity checks
    if("${PYTHON_SITE_PACKAGES_SUBDIR}" STREQUAL "")
      message(FATAL_ERROR "PYTHON_SITE_PACKAGES_SUBDIR CMake variable is expected to be set")
    endif()

    # Custom name for the components associated with ITK
    # wrapping install rules enabling Slicer to optionally
    # package ITK Wrapping in Slicer installer by simply
    # toggling the Slicer_INSTALL_ITKPython option.
    set(Slicer_WRAP_ITK_INSTALL_COMPONENT_IDENTIFIER "Wrapping")
    mark_as_superbuild(Slicer_WRAP_ITK_INSTALL_COMPONENT_IDENTIFIER:STRING)

    set(PY_SITE_PACKAGES_PATH lib/Python/${PYTHON_SITE_PACKAGES_SUBDIR})

    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}
      -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
      -DSWIG_EXECUTABLE:PATH=${SWIG_EXECUTABLE}
      -DITK_USE_SYSTEM_SWIG:BOOL=ON
      -DITK_LEGACY_SILENT:BOOL=OFF
      -DWRAP_ITK_INSTALL_COMPONENT_IDENTIFIER:STRING=${Slicer_WRAP_ITK_INSTALL_COMPONENT_IDENTIFIER}
      -DPY_SITE_PACKAGES_PATH:STRING=${PY_SITE_PACKAGES_PATH}
      )
  endif()

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  message(STATUS "Building against TBB_DIR:${TBB_DIR}:")
  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    CMAKE_CACHE_ARGS
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_EXAMPLES:BOOL=OFF
      -DITK_LEGACY_REMOVE:BOOL=ON
      -DITK_FUTURE_LEGACY_REMOVE:BOOL=ON
      -DITKV3_COMPATIBILITY:BOOL=OFF
      -DITK_BUILD_DEFAULT_MODULES:BOOL=ON
      -DModule_AnisotropicDiffusionLBR:BOOL=ON
      -DModule_ITKReview:BOOL=ON
      -DModule_MGHIO:BOOL=ON
      -DModule_ITKIOMINC:BOOL=ON
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DITK_INSTALL_NO_DEVELOPMENT:BOOL=ON
      -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
      -DITK_WRAPPING:BOOL=OFF #${BUILD_SHARED_LIBS} ## HACK:  QUICK CHANGE
      -DITK_WRAP_PYTHON:BOOL=OFF
      -DExternalData_OBJECT_STORES:PATH=${ExternalData_OBJECT_STORES}
      # macOS
      -DCMAKE_MACOSX_RPATH:BOOL=0
      # VTK
      ${ITK_VTK_OPTIONS}
      # DCMTK
      -DITK_USE_SYSTEM_DCMTK:BOOL=ON
      -DDCMTK_DIR:PATH=${DCMTK_DIR}
      -DModule_ITKIODCMTK:BOOL=ON
      # ZLIB
      -DITK_USE_SYSTEM_ZLIB:BOOL=ON
      -DZLIB_ROOT:PATH=${ZLIB_ROOT}
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
      -DITK_USE_FFTWD:BOOL=ON
      -DITK_USE_FFTWF:BOOL=ON
      -DITK_USE_GOLD_LINKER:BOOL=OFF ## RHEL7 fails to build GDCM with gold linker
      # TBB
      -DModule_ITKTBB:BOOL=ON
      -DTBB_DIR:PATH=${TBB_DIR}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(ITK_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  set(_lib_subdir lib)
  if(WIN32)
    set(_lib_subdir bin)
  endif()

  # library paths
  set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD ${ITK_DIR}/${_lib_subdir}/<CMAKE_CFG_INTDIR>)
  mark_as_superbuild(
    VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
    LABELS "LIBRARY_PATHS_LAUNCHER_BUILD"
    )

  if(Slicer_BUILD_ITKPython)
    # pythonpath
    set(${proj}_PYTHONPATH_LAUNCHER_BUILD
      ${ITK_DIR}/Wrapping/Generators/Python/<CMAKE_CFG_INTDIR>
      ${ITK_DIR}/lib/<CMAKE_CFG_INTDIR>
      ${ITK_DIR}/lib
      )
    mark_as_superbuild(
      VARS ${proj}_PYTHONPATH_LAUNCHER_BUILD
      LABELS "PYTHONPATH_LAUNCHER_BUILD"
      )
  endif()

  #-----------------------------------------------------------------------------
  # Launcher setting specific to install tree

  # Since ITK Wrapping is installed in the Slicer standard site-packages
  # location, there is no need to specify custom setting for the install
  # case.

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS ITK_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
