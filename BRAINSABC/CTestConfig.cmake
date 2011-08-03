## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   enable_testing()
##   include(CTest)
set(CTEST_PROJECT_NAME "BRAINSABC")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "testing.psychiatry.uiowa.edu")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=BRAINSABC")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_TEST_TIMEOUT 3600)  ## Set timeout to one hour for now.  There is a lot of work to be done.
