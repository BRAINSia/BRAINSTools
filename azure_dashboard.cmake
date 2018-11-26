
# set_from_env
# ------------
#
# Sets a CMake variable from an environment variable. If the
# environment variable  is not  defined then the CMake variable is not
# modified. If DEFAULT is specified then if the environment variable
# is not defined the default value is used. Alternatively, if REQUIRED
# is specified then a FATAL_ERROR is generated.
#
# set_from_env( <variable> <environment variable> [REQUIRED|DEFAULT value] )
function(set_from_env var env_var)
  if(NOT DEFINED ENV{${env_var}})
    if (ARGV2 STREQUAL "REQUIRED")
      message(FATAL_ERROR "Required environment variable \"${env_var}\" not defined.")
    elseif (ARGV2 STREQUAL "DEFAULT")
      set(${var} ${ARGV3} PARENT_SCOPE)
    endif()
  else()
    set(${var} $ENV{${env_var}} PARENT_SCOPE)
  endif()
endfunction()

set_from_env(CTEST_SITE "AGENT_MACHINENAME" REQUIRED)
set(CTEST_SITE "Azure.${CTEST_SITE}")
set(CTEST_UPDATE_VERSION_ONLY 1)

set_from_env(PARALLEL_LEVEL "PARALLEL_LEVEL" DEFAULT 3)
set(CTEST_TEST_ARGS ${CTEST_TEST_ARGS} PARALLEL_LEVEL ${PARALLEL_LEVEL})


set_from_env(workspace "AGENT_BUILDDIRECTORY" REQUIRED)
file(TO_CMAKE_PATH "${workspace}" CTEST_DASHBOARD_ROOT)

set_from_env(local_build_sourcesdirectory "BUILD_SOURCESDIRECTORY" REQUIRED)
file(RELATIVE_PATH dashboard_source_name "${workspace}" "${local_build_sourcesdirectory}")
# Make environment variables to CMake variables for CTest
set_from_env(CTEST_CMAKE_GENERATOR "CTEST_CMAKE_GENERATOR" DEFAULT "Ninja")
set_from_env(CTEST_BUILD_CONFIGURATION "CTEST_BUILD_CONFIGURATION" DEFAULT "Release")

set_from_env(BUILD_SHARED_LIBS "BUILD_SHARED_LIBS" DEFAULT "OFF")
set_from_env(BUILD_EXAMPLES "BUILD_EXAMPLES" DEFAULT "ON")

if(NOT CTEST_BUILD_NAME)
  if(DEFINED ENV{SYSTEM_PULLREQUEST_SOURCEBRANCH})
    set(branch "-$ENV{SYSTEM_PULLREQUEST_SOURCEBRANCH}")
    set(dashboard_model "Experimental")
  else()
    set(branch "-master")
  endif()

  if(DEFINED ENV{SYSTEM_PULLREQUEST_PULLREQUESTNUMBER})
    set(pr "-PR$ENV{SYSTEM_PULLREQUEST_PULLREQUESTNUMBER}")
  else()
    set(pr "")
  endif()

  set(CTEST_BUILD_NAME
    "$ENV{AGENT_OS}-Build$ENV{BUILD_BUILDID}${pr}${branch}${wrapping}")
endif()

set(dashboard_cache "
    BUILD_EXAMPLES:BOOL=${BUILD_EXAMPLES}
    BUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
BRAINSTools_BUILD_DICOM_SUPPORT:BOOL=OFF
USE_DWIConvert:BOOL=OFF
USE_ANTS:BOOL=OFF
BRAINSTools_USE_CTKAPPLAUNCHER:BOOL=OFF
BRAINSTools_USE_QT:BOOL=OFF
USE_BRAINSDWICleanup:BOOL=OFF
USE_AutoWorkup:BOOL=OFF
USE_BRAINSInitializedControlPoints:BOOL=OFF
USE_BRAINSLabelStats:BOOL=OFF
USE_BRAINSSnapShotWriter:BOOL=OFF
USE_BRAINSSuperResolution:BOOL=OFF
USE_ReferenceAtlas:BOOL=OFF
BRAINSTools_REQUIRES_VTK:BOOL=OFF
" )

string(TIMESTAMP build_date "%Y-%m-%d")
message("CDash Build Identifier: ${build_date} ${CTEST_BUILD_NAME}")
message("CTEST_SITE = ${CTEST_SITE}")
include("${CTEST_SCRIPT_DIRECTORY}/BRAINSTools_common.cmake")
