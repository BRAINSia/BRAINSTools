
#################
## Use introspection to Identify the compiler
## flags for mex compilations
#################
#
# mex -v outputs all the settings used for building MEX files, so it
# we can use it to grab the important variables needed to generate
# a well formed mex file.
execute_process(COMMAND ${MATLAB_MEX_PATH} -v ${CMAKE_CURRENT_LIST_DIR}/mex_stub.cxx
  OUTPUT_VARIABLE mexOut
  ERROR_VARIABLE mexErr)
message(STATUS "_mexOut : ${_mexOut} ")

# parse mex output.
# parse line by line by turning file into CMake list of lines
string(REGEX REPLACE "\r?\n" ";" _mexOut "${mexOut}")
message(STATUS "_mexOut : ${_mexOut} ")

foreach(line ${_mexOut})
  if("${line}" MATCHES " CXXFLAGS *:")
    string(REGEX REPLACE " *CXXFLAGS *: *" "" mexCxxFlags "${line}")
  elseif("${line}" MATCHES " CXXLIBS *:")
    string(REGEX REPLACE " *CXXLIBS *: *" "" mexCxxLibs "${line}")
  elseif("${line}" MATCHES " LDFLAGS *:")
    string(REGEX REPLACE " *LDFLAGS *: *" "" mexLdFlags "${line}")
  elseif("${line}" MATCHES " LDEXTENSION *:")
    string(REGEX REPLACE " *LDEXTENSION *: *" "" mexLdExtension "${line}")
  endif()
  message(STATUS "line:${line}")
endforeach()
list(APPEND mexCxxFlags "-DMATLAB_MEX_FILE")
message(STATUS "mexLdExtension: ${mexLdExtension}")
message(STATUS "mexLdFlags: ${mexLdFlags}")
message(STATUS "mexCxxLibs: ${mexCxxLibs}")
message(FATAL_ERROR "mexCxxFlags: ${mexCxxFlags}")

#####################
#####################
/* The following is the world's smallest MEX-file that
   actually does something.
   Note that it is written in ANSI C.
   */
#include "mex.h"

void mexFunction(int nlhs,
  Matrix *plhs[],
  int nrhs,
  Matrix *prhs[])
{
  mexPrintf("Hello world");
}
