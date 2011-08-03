#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
# #-----------------------------------------------------------------------------
if(POLICY CMP0016)
  cmake_policy(SET CMP0016 NEW)
endif()
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

include(ExternalProject)

if(INTEGRATE_WITH_SLICER)
  set(BRAINS_BUILD OFF CACHE INTERNAL "Set BRAINS_BUILD to off for slicer builds" FORCE )
else()
  set(BRAINS_BUILD ON CACHE INTERNAL "Set BRAINS_BUILD to on for non-slicer builds" FORCE )
endif()

#

#-----------------------------------------------------------------------------
# Build the optional DEBUGIMAGEVIEWER
if(NOT SETOPTIONALDEBUGIMAGEVIEWER)
macro(SETOPTIONALDEBUGIMAGEVIEWER)
if(BRAINS_BUILD)
  option(USE_DEBUG_IMAGE_VIEWER "Use the DEBUG_IMAGE_VIEWER for debugging" ON)
else(BRAINS_BUILD)
  option(USE_DEBUG_IMAGE_VIEWER "Use the DEBUG_IMAGE_VIEWER for debugging" OFF)
endif(BRAINS_BUILD)

mark_as_advanced(USE_DEBUG_IMAGE_VIEWER)
set(OPTIONAL_DEBUG_LINK_LIBRARIES) ## Set it to empty as the default
if( USE_DEBUG_IMAGE_VIEWER )
   if(NOT KWWidgets_SOURCE_DIR)
     find_package(KWWidgets REQUIRED)
     include(${KWWidgets_USE_FILE})
   endif(NOT KWWidgets_SOURCE_DIR)
   add_definitions(-DUSE_DEBUG_IMAGE_VIEWER)
   find_path(DEBUG_IMAGE_VIEWER_INCLUDE_DIR DebugImageViewerClient.h ${CMAKE_INSTALL_PREFIX}/include)
   include_directories(${DEBUG_IMAGE_VIEWER_INCLUDE_DIR})
   set(OPTIONAL_DEBUG_LINK_LIBRARIES ${KWWidgets_LIBRARIES})
endif( USE_DEBUG_IMAGE_VIEWER )
endmacro(SETOPTIONALDEBUGIMAGEVIEWER)
endif(NOT SETOPTIONALDEBUGIMAGEVIEWER)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## A macro to create executables for Slicer or BRAINS3
if(NOT CONFIGUREBRAINSORSLICERPROPERTIES)
  macro(CONFIGUREBRAINSORSLICERPROPERTIES PROGNAME PROGCLI PROGSOURCES LIBSOURCES ENTRYPOINTNAME EXTRA_LIBS)

  find_package(SlicerExecutionModel NO_MODULE REQUIRED GenerateCLP)
  include(${GenerateCLP_USE_FILE})

  get_filename_component(TMP_FILENAME ${PROGCLI} NAME_WE)
  set(PROGCLI_HEADER "${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h")

  set(CLP_SOURCES ${PROGSOURCES} ${LIBSOURCES})
  set(CLP_PRIMARY_SOURCES ${LIBSOURCES})
  if(EXISTS  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h)
    GENERATECLP(CLP_SOURCES ${PROGCLI} ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h)
  else()
    GENERATECLP(CLP_SOURCES ${PROGCLI} )
  endif()

  add_executable( ${PROGNAME} ${CLP_SOURCES} ${PROGCLI_HEADER})

  if(WIN32)
    set(BRAINS_ITK_LIBS "")
  else(WIN32)
    if(ITK_VERSION_MAJOR LESS "4")
      set(BRAINS_ITK_LIBS ITKAlgorithms ITKIO ITKBasicFilters)
    else()
      set(BRAINS_ITK_LIBS ${ITK_LIBRARIES})
    endif()
  endif(WIN32)
  target_link_libraries (${PROGNAME} BRAINSCommonLib ${BRAINS_ITK_LIBS} ${OPTIONAL_DEBUG_LINK_LIBRARIES} ${EXTRA_LIBS} )

  if (INTEGRATE_WITH_SLICER)
    #if(NOT slicer3_set_plugins_output_path)
    #  include(${Slicer_SOURCE_DIR}/CMake/Slicer3PluginsMacros.cmake)
    #endif()
    ### If building as part of the INTEGRATE_WITH_SLICER, then only build the shared object, and not the command line program.

    message(STATUS "Building ${PROGNAME} for Slicer4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    add_library(${PROGNAME}Lib SHARED ${CLP_PRIMARY_SOURCES} ${PROGCLI_HEADER})
    set_target_properties (${PROGNAME}Lib PROPERTIES COMPILE_FLAGS "-D${ENTRYPOINTNAME}=ModuleEntryPoint")
    target_link_libraries (${PROGNAME}Lib BRAINSCommonLib ${BRAINS_ITK_LIBS} ${OPTIONAL_DEBUG_LINK_LIBRARIES} ${EXTRA_LIBS} )

    # install each target in the production area (where it would appear in an
    # installation) and install each target in the developer area (for running
    # from a build)
    set(TARGETS ${PROGNAME}Lib ${PROGNAME})
    #slicer3_set_plugins_output_path(${PROGNAME}Lib)
    #slicer3_set_plugins_output_path(${PROGNAME})
    #slicer3_install_plugins(${TARGETS})
    install(TARGETS ${PROGNAME}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            COMPONENT Development
            )

  else (INTEGRATE_WITH_SLICER)
  endif (INTEGRATE_WITH_SLICER)
    ### If building outside of Slicer3, then only build the command line executable.
    install(TARGETS ${PROGNAME}
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      COMPONENT Development)

endmacro(CONFIGUREBRAINSORSLICERPROPERTIES PROGNAME PROGCLI PROGSOURCES LIBSOURCES)
endif(NOT CONFIGUREBRAINSORSLICERPROPERTIES)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## A macro to create CLP dependant libraries for Slicer or BRAINS3
if(NOT CONFIGUREBRAINSORSLICERLIBRARY)
macro(CONFIGUREBRAINSORSLICERLIBRARY LIBNAME LIBCLI LIBSOURCES EXTRA_LIBS)

  set(CLP_SOURCES ${LIBSOURCES})
  if("x${LIBCLI}" MATCHES "x")
    ## Nothing to do if there is no CLI component
  else()
    find_package(SlicerExecutionModel NO_MODULE REQUIRED GenerateCLP)
    include(${GenerateCLP_USE_FILE})

    get_filename_component(TMP_FILENAME ${LIBCLI} NAME_WE)
    set(LIBCLI_HEADER "${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h")

    if(EXISTS  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h)
      GENERATECLP(CLP_SOURCES ${LIBCLI} ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h)
    else()
      GENERATECLP(CLP_SOURCES ${LIBCLI} )
    endif()
  endif()

  add_library( ${LIBNAME} ${CLP_SOURCES} ${LIBCLI_HEADER})
  target_link_libraries (${LIBNAME} BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES} ${EXTRA_LIBS} )

  ##HACK:  if (INTEGRATE_WITH_SLICER)
    #if(NOT slicer3_set_plugins_output_path)
    #  include(${Slicer_SOURCE_DIR}/CMake/Slicer3PluginsMacros.cmake)
    #endif()
    ### If building as part of the INTEGRATE_WITH_SLICER, then only build the shared object, and not the command line program.
    # install each target in the production area (where it would appear in an
    # installation) and install each target in the developer area (for running
    # from a build)
    #slicer3_set_plugins_output_path(${LIBNAME})
    #set(TARGETS ${LIBNAME})
    #slicer3_install_plugins(${TARGETS})
  ##HACK:  else (INTEGRATE_WITH_SLICER)
    install(TARGETS ${LIBNAME}
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      COMPONENT Development)
  ##HACK:  endif (INTEGRATE_WITH_SLICER)

endmacro(CONFIGUREBRAINSORSLICERLIBRARY LIBNAME LIBCLI LIBSOURCES)
endif(NOT CONFIGUREBRAINSORSLICERLIBRARY)
