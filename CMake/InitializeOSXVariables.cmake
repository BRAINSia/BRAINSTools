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

  # Disable universal binaries
  if(NOT "${CMAKE_OSX_ARCHITECTURES}" STREQUAL "")
    list(LENGTH CMAKE_OSX_ARCHITECTURES arch_count)
    if(arch_count GREATER 1)
      message(FATAL_ERROR "error: Only one value (arm64, i386 or x86_64) should be associated with CMAKE_OSX_ARCHITECTURES.")
    endif()
  endif()

  set(minimum_allowed_deployment_target "12.0")
  if("x${CMAKE_OSX_DEPLOYMENT_TARGET}x" STREQUAL "xx")
    execute_process(
      COMMAND xcrun --sdk macosx --show-sdk-version
      OUTPUT_VARIABLE autodected_deployment_target
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Detected deployment target: ${autodected_deployment_target}")
    set(CMAKE_OSX_DEPLOYMENT_TARGET ${autodected_deployment_target} CACHE STRING "OSX Deployment target" FORCE)

    # Get the SDK path using xcrun
    execute_process(
      COMMAND xcrun --sdk macosx --show-sdk-path
      OUTPUT_VARIABLE CMAKE_OSX_SYSROOT
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Optionally set it in the cache for consistency
    set(CMAKE_OSX_SYSROOT "${CMAKE_OSX_SYSROOT}" CACHE PATH "Path to the macOS SDK")
    message(STATUS "Detected macOS SDK: ${CMAKE_OSX_SYSROOT}")
  endif()

  if(CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS ${minimum_allowed_deployment_target})
    message(FATAL_ERROR "CMAKE_OSX_DEPLOYMENT_TARGET ${CMAKE_OSX_DEPLOYMENT_TARGET} must be ${minimum_allowed_deployment_target} or greater.")
  endif()

  if(NOT "${CMAKE_OSX_SYSROOT}" STREQUAL "")
    if(NOT EXISTS "${CMAKE_OSX_SYSROOT}")
      message(FATAL_ERROR "error: CMAKE_OSX_SYSROOT='${CMAKE_OSX_SYSROOT}' does not exist")
    endif()
  endif()
endif()
