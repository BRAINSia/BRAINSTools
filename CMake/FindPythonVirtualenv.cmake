#.rst
# FindPythonVirtualenv
# --------------------
#
# Find python and libraries in a virtualenv
#
# Virtualenv should be built with the --always-copy and --python=${PYTHONEXEC} flags set
# and the python-config file should be copied to the virtualenv bin/ directory
#
# This module determines if build is being done within a virtualenv and where the Python
# executable, include files, and libraries are according to the $VIRTUAL_ENV. Calls
# find_package( PythonInterp) and if no virtualenv, calls find_package( PythonLibs). This
# code sets the following variables:
#
# ::
#
#   PYTHON_EXECUTABLE          - path to the virtualenv Python interpreter
#   PYTHON_LIBRARY             - path to the virtualenv Python library
#   PYTHON_INCLUDE_DIR         - path to where virtualenv Python.h is found
#   PYTHON_VERSION_STRING      - Virtualenv Python version found, e.g. 2.5.2
#   PYTHON_VERSION_MAJOR       - Python major version found, e.g. 2
#   PYTHON_VERSION_MINOR       - Python minor version found, e.g. 5
#   PYTHON_VERSION_PATCH       - Python patch version found, e.g. 2
#
#


find_package( PythonInterp REQUIRED )

if( DEFINED ENV{VIRTUAL_ENV} )
  message(STATUS "Found Virtualenv: $ENV{VIRTUAL_ENV}")
  # VirtualEnv MUST be activated
  if( IS_DIRECTORY "$ENV{VIRTUAL_ENV}" )
    # Overwrite the executable found in FindPythonInterp.cmake, if needed
    set( PYTHON_EXECUTABLE "$ENV{VIRTUAL_ENV}/bin/python" CACHE FILEPATH "Path to Python virtualenv executable" FORCE)
    # Update python version string
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                            "import sys; sys.stdout.write(';'.join([str(x) for x in sys.version_info[:3]]))"
                    OUTPUT_VARIABLE _VERSION
                    RESULT_VARIABLE _PYTHON_VERSION_RESULT
                    ERROR_QUIET)
    if(NOT _PYTHON_VERSION_RESULT)
      string(REPLACE ";" "." PYTHON_VERSION_STRING "${_VERSION}")
      list(GET _VERSION 0 PYTHON_VERSION_MAJOR)
      list(GET _VERSION 1 PYTHON_VERSION_MINOR)
      list(GET _VERSION 2 PYTHON_VERSION_PATCH)
      if(PYTHON_VERSION_PATCH EQUAL 0)
        # it's called "Python 2.7", not "2.7.0"
        string(REGEX REPLACE "\\.0$" "" PYTHON_VERSION_STRING "${PYTHON_VERSION_STRING}")
      endif()
    else()
      # sys.version predates sys.version_info, so use that
      execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import sys; sys.stdout.write(sys.version)"
        OUTPUT_VARIABLE _VERSION
        RESULT_VARIABLE _PYTHON_VERSION_RESULT
        ERROR_QUIET)
      if(NOT _PYTHON_VERSION_RESULT)
        string(REGEX REPLACE " .*" "" PYTHON_VERSION_STRING "${_VERSION}")
        string(REGEX REPLACE "^([0-9]+)\\.[0-9]+.*" "\\1" PYTHON_VERSION_MAJOR "${PYTHON_VERSION_STRING}")
        string(REGEX REPLACE "^[0-9]+\\.([0-9])+.*" "\\1" PYTHON_VERSION_MINOR "${PYTHON_VERSION_STRING}")
        if(PYTHON_VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.([0-9]+)")
          set(PYTHON_VERSION_PATCH "${CMAKE_MATCH_1}")
        else()
          set(PYTHON_VERSION_PATCH "0")
        endif()
      else()
        # sys.version was first documented for Python 1.5, so assume
        # this is older.
        set(PYTHON_VERSION_STRING "1.4")
        set(PYTHON_VERSION_MAJOR "1")
        set(PYTHON_VERSION_MINOR "4")
        set(PYTHON_VERSION_PATCH "0")
      endif()
    endif()
    unset(_PYTHON_VERSION_RESULT)
    unset(_VERSION)
  endif()
  # Set values for FindPythonLibs.cmake
  string(REGEX REPLACE "\\.[0-9]$" "" PYTHON_MAJOR_MINOR "${PYTHON_VERSION_STRING}")
  set( PYTHON_LIBRARY "$ENV{VIRTUAL_ENV}/lib/python${PYTHON_MAJOR_MINOR}" CACHE PATH "Path to Python virtualenv libs" FORCE )
  # message( STATUS "This is the current library: ${PYTHON_LIBRARY}" )
  set( PYTHON_INCLUDE_DIR "$ENV{VIRTUAL_ENV}/include/" CACHE PATH "Path to Python virtualenv include" FORCE )
  # message( STATUS "This is the current headers: ${PYTHON_INCLUDE_DIR}" )
  unset( PYTHON_MAJOR_MINOR )
endif()

find_package( PythonLibs REQUIRED )

mark_as_advanced(
  PYTHON_EXECUTABLE
  PYTHON_LIBRARY
  PYTHON_INCLUDE_DIR
  )

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  PythonVirtualenv
  REQUIRED_VARS PYTHON_EXECUTABLE PYTHON_LIBRARY PYTHON_INCLUDE_DIR
  VERSION_VAR PYTHON_VERSION_STRING
  FAIL_MESSAGE "FindPythonVirtualenv failed to find complete Python environment")
