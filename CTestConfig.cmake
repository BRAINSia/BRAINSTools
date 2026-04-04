## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   enable_testing()
##   include(CTest)

set(CTEST_PROJECT_NAME "BRAINSTools")
set(CTEST_NIGHTLY_START_TIME "00:00:00 UTC")

# CTEST_SUBMIT_URL (CMake >= 3.14) consolidates the legacy CTEST_DROP_*
# variables into one setting.  GitHub Actions runners ship CMake >= 3.25
# so the legacy fallback is retained only for local builds with older CMake.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.14")
  set(CTEST_SUBMIT_URL "https://my.cdash.org/submit.php?project=BRAINSTools")
else()
  set(CTEST_DROP_METHOD "https")
  set(CTEST_DROP_SITE "my.cdash.org")
  set(CTEST_DROP_LOCATION "/submit.php?project=BRAINSTools")
endif()

set(CTEST_DROP_SITE_CDASH TRUE)
