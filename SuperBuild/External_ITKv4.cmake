
set(proj ITKv4)
set(ITK_EXTERNAL_NAME ${proj})
# Set dependency list
set(${proj}_DEPENDENCIES "zlib" VTK)
#if(${CMAKE_PROJECT_NAME}_BUILD_DICOM_SUPPORT)
  list(APPEND ${proj}_DEPENDENCIES DCMTK)
#endif()

if(${PRIMARY_PROJECT_NAME}_USE_QT) ## QT requires VTK support in ITK
  list(APPEND ${proj}_DEPENDENCIES VTK)
endif()

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(ITK_DIR CACHE)
  find_package(ITK 4 COMPONENTS ${${CMAKE_PROJECT_NAME}_ITK_COMPONENTS} REQUIRED NO_MODULE)
endif()

# Sanity checks
if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to non-existing directory")
endif()

if(NOT DEFINED ITK_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  if(BRAINSTools_REQUIRES_VTK) ## QT requires VTK support in ITK
    set( ITK_VTK_OPTIONS
      -DVTK_DIR:PATH=${VTK_DIR}
      -DModule_ITKVtkGlue:BOOL=ON  ## If building with GUI, then need ITKVtkGlue
    )
  endif()

  if(NOT DEFINED git_protocol)
      set(git_protocol "git")
  endif()

  #set(${proj}_REPOSITORY ${git_protocol}://itk.org/ITK.git)
  set(${proj}_REPOSITORY ${git_protocol}://github.com/BRAINSia/ITK.git)
  set(${proj}_GIT_TAG c123bf502918db7acb7b9afe5ba3fc5cd362b47e )
#  message("COMMON_EXTERNAL_PROJECT_ARGS:
#${COMMON_EXTERNAL_PROJECT_ARGS}")
  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${ITKv4_REPOSITORY}
    GIT_TAG ${ITKv4_GIT_TAG}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${proj}-build
    CMAKE_CACHE_ARGS
      ${COMMON_EXTERNAL_PROJECT_ARGS}
    CMAKE_ARGS
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_EXAMPLES:BOOL=OFF
      -DITK_LEGACY_REMOVE:BOOL=OFF
      -DITKV3_COMPATIBILITY:BOOL=ON
      -DITK_BUILD_DEFAULT_MODULES:BOOL=ON
      -DModule_ITKReview:BOOL=ON
      ${ITK_VTK_OPTIONS}
      -DModule_MGHIO:BOOL=ON        #To provide FreeSurfer Compatibility
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DITK_INSTALL_NO_DEVELOPMENT:BOOL=ON
      -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
      -DITK_WRAPPING:BOOL=OFF #${BUILD_SHARED_LIBS} ## HACK:  QUICK CHANGE
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
      -DVTK_DIR:PATH=${VTK_DIR}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )
 set(ITK_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS ITK_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
