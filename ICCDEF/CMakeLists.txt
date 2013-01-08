##- project(iccdefRegistrationNew)

#-------------------------------------------------------------------------------
# OUTPUT DIRECTORIES
#
##- include_directories(${COMMON_BRAINSCOMMONLIB_SOURCE_DIR})
##- link_directories(${COMMON_BRAINSCOMMONLIB_BINARY_DIR})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/iccdefRegistration_New" CACHE PATH "Install path prefix,
     prepended onto install directories." FORCE)
endif()

mark_as_advanced(ITK_USE_FFTWD)
option(ITK_USE_FFTWD "Use double precision FFTW if found" OFF)
mark_as_advanced(ITK_USE_FFTWF)
option(ITK_USE_FFTWF "Use single precision FFTW if found" OFF)

if( ITK_USE_FFTWD OR ITK_USE_FFTWF )
    find_package( FFTW )
    if( ITK_USE_FFTWD )
      if( NOT FFTWD_FOUND)
        message(FATAL_ERROR "Cannot build iccdef without FFTW. Please set FFTWD_INCLUDE and FFTWD_LIB")
      endif()
        link_libraries(${FFTWD_LIB})
     endif()

    if( ITK_USE_FFTWF )
      if( NOT FFTWF_FOUND )
        message(FATAL_ERROR "Cannot build iccdef without FFTW. Please set FFTWF_INCLUDE and FFTWF_LIB")
      endif()
        link_libraries(${FFTWF_LIB})
    endif()
endif()

include_directories(
      ${FFTW_INCLUDE}
)

set(ALL_PROGS_LIST
  ApplyWarp
  AverageBrainGenerator
  CreateMask
  iccdefRegistration_New
  )

set(ICCDEFLibraries BRAINSCommonLib  ${FFTWF_LIB})
DebugImageViewerLibAdditions(ICCDEFLibraries)

foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES ${ICCDEFLibraries})
endforeach()

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
