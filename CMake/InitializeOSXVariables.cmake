################################################################################
#
#  Program: 3D Slicer
#
#  Copyright (c) Kitware Inc.
#
#  See COPYRIGHT.txt
#  or http://www.slicer.org/copyright/copyright.txt for details.
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
#  and was partially funded by NIH grant 3P41RR013218-12S1
#
################################################################################

#
# SlicerInitializeOSXVariables
#

#
# Adapted from Paraview/Superbuild/CMakeLists.txt
#

# Note: Change architecture *before* any enable_language() or project()
#       calls so that it's set properly to detect 64-bit-ness, and
#       deployment target for the standard c++ library.
#
if(APPLE)
  #----------------------------------------------------------------------------
  # _CURRENT_OSX_VERSION - as a two-component string: 10.11, 10.12, 10.13, 10.14 ...
  #
  string(REGEX REPLACE "^([0-9]+\\.[0-9]+).*$" "\\1"
    _CURRENT_OSX_VERSION "${CURRENT_OSX_VERSION}")

  # Waiting universal binaries are supported and tested, complain if
  # multiple architectures are specified.
  if(NOT "${CMAKE_OSX_ARCHITECTURES}" STREQUAL "")
    list(LENGTH CMAKE_OSX_ARCHITECTURES arch_count)
    if(arch_count GREATER 1)
      message(FATAL_ERROR "error: Only one value (i386 or x86_64) should be associated with CMAKE_OSX_ARCHITECTURES.")
    endif()
  endif()

  # See CMake/Modules/Platform/Darwin.cmake and https://en.wikipedia.org/wiki/MacOS#Release_history
  #   8.x == Mac OSX 10.4 (Tiger)
  #   9.x == Mac OSX 10.5 (Leopard)
  #  10.x == Mac OSX 10.6 (Snow Leopard)
  #  11.x == Mac OSX 10.7 (Lion)
  #  12.x == Mac OSX 10.8 (Mountain Lion)
  #  13.x == Mac OSX 10.9 (Mavericks)
  #  14.x == Mac OSX 10.10 (Yosemite)
  #  15.x == Mac OSX 10.11 (El Capitan)
  #  16.x == Mac OSX 10.12 (Sierra)
  #  17.x == Mac OSX 10.13 (High Sierra)
  #  18.x == Mac OSX 10.14 (Mojave)
  set(OSX_SDK_104_NAME "Tiger")
  set(OSX_SDK_105_NAME "Leopard")
  set(OSX_SDK_106_NAME "Snow Leopard")
  set(OSX_SDK_107_NAME "Lion")
  set(OSX_SDK_108_NAME "Mountain Lion")
  set(OSX_SDK_109_NAME "Mavericks")
  set(OSX_SDK_1010_NAME "Yosemite")
  set(OSX_SDK_1011_NAME "El Capitan")
  set(OSX_SDK_1012_NAME "Sierra")
  set(OSX_SDK_1013_NAME "High Sierra")
  set(OSX_SDK_1014_NAME "Mojave")

  set(OSX_SDK_ROOTS
    /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs
    /Developer/SDKs
    )

  # Explicitly set the OSX_SYSROOT to the latest one, as its required
  #       when the SX_DEPLOYMENT_TARGET is explicitly set
  foreach(SDK_ROOT ${OSX_SDK_ROOTS})
    if( "x${CMAKE_OSX_SYSROOT}x" STREQUAL "xx")
      file(GLOB SDK_SYSROOTS "${SDK_ROOT}/MacOSX*.sdk")

      if(NOT "x${SDK_SYSROOTS}x" STREQUAL "xx")
        set(SDK_SYSROOT_NEWEST "")
        set(SDK_VERSION "0")
        # find the latest SDK
        foreach(SDK_ROOT_I ${SDK_SYSROOTS})
          # extract version from SDK
          string(REGEX MATCH "MacOSX([0-9]+\\.[0-9]+)\\.sdk" _match "${SDK_ROOT_I}")
          if("${CMAKE_MATCH_1}" VERSION_GREATER "${SDK_VERSION}")
            set(SDK_SYSROOT_NEWEST "${SDK_ROOT_I}")
            set(SDK_VERSION "${CMAKE_MATCH_1}")
          endif()
        endforeach()

        if(NOT "x${SDK_SYSROOT_NEWEST}x" STREQUAL "xx")
          string(REPLACE "." "" sdk_version_no_dot ${SDK_VERSION})
          set(OSX_NAME ${OSX_SDK_${sdk_version_no_dot}_NAME})
          set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Force build for 64-bit ${OSX_NAME}." FORCE)
          set(CMAKE_OSX_SYSROOT "${SDK_SYSROOT_NEWEST}" CACHE PATH "Force build for 64-bit ${OSX_NAME}." FORCE)
          message(STATUS "Setting OSX_ARCHITECTURES to '${CMAKE_OSX_ARCHITECTURES}' as none was specified.")
          message(STATUS "Setting OSX_SYSROOT to latest '${CMAKE_OSX_SYSROOT}' as none was specified.")
        endif()
      endif()
    endif()
  endforeach()

  if(NOT "${CMAKE_OSX_SYSROOT}" STREQUAL "")
    if(NOT EXISTS "${CMAKE_OSX_SYSROOT}")
      message(FATAL_ERROR "error: CMAKE_OSX_SYSROOT='${CMAKE_OSX_SYSROOT}' does not exist")
    endif()
  endif()

    # In 10.9 libc++ replaces libstdc++
    # as the default runtime. Requiring this minimum ensures that all libraries
    # use libc++.
    set(required_deployment_target "10.13")
  ## See https://doc.qt.io/qt-5.11/supported-platforms-and-configurations.html for supported versions of qt
  ## Provide backwards compatibility levels equal to QT support (i.e. 10.14 QT 5.12 also supports 10.13, and 10.13)
  if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    if( _CURRENT_OSX_VERSION VERSION_GREATER "10.13" )
       set(CMAKE_OSX_DEPLOYMENT_TARGET "${_CURRENT_OSX_VERSION}" CACHE STRING "Force deployment target to 10.14" FORCE)
    else()
       set(CMAKE_OSX_DEPLOYMENT_TARGET "${required_deployment_target}" CACHE STRING "Force deployment target to 10.14" FORCE)
    endif()
  endif()

  if(CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS ${required_deployment_target})
    message(FATAL_ERROR "CMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET} must be ${required_deployment_target} or greater.")
  endif()
  if("x${CMAKE_OSX_DEPLOYMENT_TARGET}x" STREQUAL "xx")
    string(REGEX MATCH "MacOSX([0-9]+\\.[0-9]+)\\.sdk" _match "${CMAKE_OSX_SYSROOT}")
    set(SDK_VERSION "${CMAKE_MATCH_1}")
    if( "${SDK_VERSION}" VERSION_LESS "10.13" )
      message(FATAL_ERROR
        "The ${SDK_VERSION} (<10.13) and OSX_DEPLOYMENT_TARGET are not supported!\n"
        )
    endif()
  endif()
endif()
