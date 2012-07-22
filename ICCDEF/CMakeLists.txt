project(iccdefRegistrationNew)
set(LOCAL_RROJECT_NAME iccdefRegistrationNew)

#-------------------------------------------------------------------------------
# OUTPUT DIRECTORIES
#
include_directories(${COMMON_BRAINSCOMMONLIB_SOURCE_DIR})
link_directories(${COMMON_BRAINSCOMMONLIB_BINARY_DIR})

set(CMAKE_MODULE_PATH ${iccdefRegistrationNew_SOURCE_DIR} ${CMAKE_MODULE_PATH} )
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/iccdefRegistration_New" CACHE PATH "Install path prefix,
     prepended onto install directories." FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

mark_as_advanced(USE_FFTWD)
option(USE_FFTWD "Use double precision FFTW if found" OFF)
mark_as_advanced(USE_FFTWF)
option(USE_FFTWF "Use single precision FFTW if found" OFF)

if( USE_FFTWD OR USE_FFTWF )
    include(${iccdefRegistrationNew_SOURCE_DIR}/FindFFTW.cmake)
    find_package( FFTW )
    if( USE_FFTWD )
      if( NOT FFTWD_FOUND)
        message(FATAL_ERROR "Cannot build iccdef without FFTW. Please set FFTWD_INCLUDE and FFTWD_LIB")
      endif(NOT FFTWD_FOUND)
        link_libraries(${FFTWD_LIB})
     endif( USE_FFTWD )

    if( USE_FFTWF )
      if( NOT FFTWF_FOUND )
        message(FATAL_ERROR "Cannot build iccdef without FFTW. Please set FFTWF_INCLUDE and FFTWF_LIB")
      endif( NOT FFTWF_FOUND )
        link_libraries(${FFTWF_LIB})
    endif( USE_FFTWF )
endif( USE_FFTWD OR USE_FFTWF )

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

if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif(BUILD_TESTING)
