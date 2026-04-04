#-----------------------------------------------------------------------------
# BRAINSToolsInstrumentation.cmake
#
# Performance monitoring and runtime-instrumentation compiler flags.
#
# Provides four mutually-exclusive options:
#   BUILD_PROFILING  – frame-pointer + debug info for Instruments/perf/gprof
#   BUILD_ASAN       – AddressSanitizer (memory errors)
#   BUILD_UBSAN      – UndefinedBehaviourSanitizer (UB at runtime)
#   BUILD_TSAN       – ThreadSanitizer (data races)
#
# Only one of these (or BUILD_COVERAGE / BUILD_OPTIMIZED) may be ON at once.
# Flags are applied only to BRAINSTools itself; ExternalProject dependencies
# (ITK, VTK, etc.) are intentionally left uninstrumented.
#
# Usage:
#   include(BRAINSToolsInstrumentation)   # from Common.cmake
#
# See also:
#   Utilities/misc/bt-configure-inner.sh  --cmake-args "-DBUILD_PROFILING=ON"
#   Utilities/misc/bt-profile-tests.sh
#-----------------------------------------------------------------------------

option(BUILD_PROFILING
  "Build with -fno-omit-frame-pointer and debug symbols for use with Instruments/perf/gprof."
  OFF)
option(BUILD_ASAN
  "Build with AddressSanitizer (-fsanitize=address). Requires clang or gcc >= 4.8."
  OFF)
option(BUILD_UBSAN
  "Build with UndefinedBehaviourSanitizer (-fsanitize=undefined)."
  OFF)
option(BUILD_TSAN
  "Build with ThreadSanitizer (-fsanitize=thread). Cannot be combined with ASAN."
  OFF)
mark_as_advanced(BUILD_PROFILING BUILD_ASAN BUILD_UBSAN BUILD_TSAN)

# Mutual-exclusion enforcement
set(_profiling_options_enabled 0)
foreach(_opt BUILD_PROFILING BUILD_ASAN BUILD_UBSAN BUILD_TSAN BUILD_COVERAGE BUILD_OPTIMIZED)
  if(${_opt})
    math(EXPR _profiling_options_enabled "${_profiling_options_enabled} + 1")
  endif()
endforeach()
if(_profiling_options_enabled GREATER 1)
  message(FATAL_ERROR
    "At most one of BUILD_PROFILING / BUILD_ASAN / BUILD_UBSAN / BUILD_TSAN / "
    "BUILD_COVERAGE / BUILD_OPTIMIZED may be ON at a time.")
endif()
unset(_profiling_options_enabled)

if(BUILD_PROFILING)
  # -fno-omit-frame-pointer: preserve call-chain frames visible to profilers.
  # -fno-optimize-sibling-calls: prevent tail-call elimination hiding callers.
  # -g: include DWARF debug info so Instruments/perf shows symbol names.
  # These are deliberately NOT forwarded to ExternalProject deps (ITK, VTK, etc.)
  # because profiling overhead in deps is rarely useful and slows EP builds.
  foreach(_lang C CXX)
    string(APPEND CMAKE_${_lang}_FLAGS
      " -fno-omit-frame-pointer -fno-optimize-sibling-calls -g")
  endforeach()
  message(STATUS
    "BUILD_PROFILING: added -fno-omit-frame-pointer -fno-optimize-sibling-calls -g")
endif()

if(BUILD_ASAN)
  # AddressSanitizer: detect heap/stack/global buffer overflows, use-after-free.
  # Requires a matching runtime library; on macOS: clang ships it in Xcode.
  # NOTE: ASAN is intentionally NOT forwarded to EP deps — instrumented ITK/VTK
  # would produce false positives from their own allocators.
  foreach(_lang C CXX)
    string(APPEND CMAKE_${_lang}_FLAGS " -fsanitize=address -fno-omit-frame-pointer -g")
  endforeach()
  string(APPEND CMAKE_EXE_LINKER_FLAGS    " -fsanitize=address")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -fsanitize=address")
  message(STATUS "BUILD_ASAN: added -fsanitize=address")
endif()

if(BUILD_UBSAN)
  foreach(_lang C CXX)
    string(APPEND CMAKE_${_lang}_FLAGS " -fsanitize=undefined -fno-omit-frame-pointer -g")
  endforeach()
  string(APPEND CMAKE_EXE_LINKER_FLAGS    " -fsanitize=undefined")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -fsanitize=undefined")
  message(STATUS "BUILD_UBSAN: added -fsanitize=undefined")
endif()

if(BUILD_TSAN)
  if(BUILD_ASAN)
    message(FATAL_ERROR "BUILD_TSAN and BUILD_ASAN cannot both be ON (incompatible runtimes).")
  endif()
  foreach(_lang C CXX)
    string(APPEND CMAKE_${_lang}_FLAGS " -fsanitize=thread -fno-omit-frame-pointer -g")
  endforeach()
  string(APPEND CMAKE_EXE_LINKER_FLAGS    " -fsanitize=thread")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -fsanitize=thread")
  message(STATUS "BUILD_TSAN: added -fsanitize=thread")
endif()
