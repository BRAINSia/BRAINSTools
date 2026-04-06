#
# BRAINSToolsGTest.cmake
# ----------------------
# Provides the macro  brainstools_enable_gtest()  which makes the CMake
# imported targets  GTest::gtest  and  GTest::gtest_main  available to
# any BRAINSTools sub-module that needs to build GTest-based unit tests.
#
# Strategy (least invasive — no SuperBuild change required):
#   1. Try find_package(GTest) to pick up a system or Homebrew installation.
#   2. If not found, fall back to FetchContent to download GTest v1.17.0
#      from GitHub.  This matches the version bundled in ITK's ThirdParty
#      module (Modules/ThirdParty/GoogleTest) so behaviour is consistent.
#
# After calling this macro, downstream CMakeLists.txt can:
#
#   add_executable(MyTest MyTest.cxx)
#   target_link_libraries(MyTest PRIVATE GTest::gtest_main)
#   include(GoogleTest)
#   gtest_discover_tests(MyTest DISCOVERY_TIMEOUT 120)
#
# The macro is guarded so it is safe to call from multiple sub-modules;
# only the first call performs the find/fetch work.
#
# Requirements:
#   CMake >= 3.20.6  (BRAINSTools minimum; FetchContent and gtest_discover_tests
#                     are fully supported from 3.14 and 3.10 respectively)
#
cmake_minimum_required(VERSION 3.20.6)

macro(brainstools_enable_gtest)
  # Guard: if the targets already exist (called by another sub-module first,
  # or the user has a find_package(GTest) higher up) do nothing.
  if(NOT TARGET GTest::gtest)
    # --- pass 1: honour an existing system/Homebrew installation ----------
    find_package(GTest QUIET)

    if(NOT GTest_FOUND)
      # --- pass 2: fetch GTest from GitHub --------------------------------
      # Version 1.17.0 matches what ITK 5.4 bundles internally.
      # GIT_SHALLOW keeps the download minimal (no history).
      # INSTALL_GTEST=OFF prevents GTest from polluting the install tree.
      # BUILD_GMOCK=OFF avoids compiling the mocking library we don't need.
      include(FetchContent)
      FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.17.0 # matches ITK's bundled GoogleTest version
        GIT_SHALLOW ON
      )
      set(INSTALL_GTEST
          OFF
          CACHE BOOL "Do not install GTest alongside BRAINSTools" FORCE)
      set(BUILD_GMOCK
          OFF
          CACHE BOOL "Do not build GMock" FORCE)
      # Suppress GTest's own install() rules so they don't pollute our
      # install tree.  FetchContent does NOT call install() for sub-projects
      # by default when EXCLUDE_FROM_ALL is set, but setting INSTALL_GTEST=OFF
      # is the official GTest-recommended guard.
      FetchContent_MakeAvailable(googletest)
    endif()
  endif()
endmacro()
