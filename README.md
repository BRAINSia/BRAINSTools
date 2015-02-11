The BRAINSTools is a harness to assist in building the many of the BRAINSTools under development.

Developers should run the "./Utilities/SetupForDevelopment.sh" script to get started.

If developing on Mac OS X, make sure the xcode command line tools are installed with the command:
xcode-select --install

For more information on the individual BRAINSTools, please see the following link:
https://github.com/BRAINSia/BRAINSTools/wiki

## Building ##
Example session for a clean build:
mkdir ~/src/mytester
cd ~/src/mytester
git clone git://github.com/BRAINSia/BRAINSTools.git
cd BRAINSTools/
bash ./Utilities/SetupForDevelopment.sh
mkdir -p ../BRAINSTools-build
cd ../BRAINSTools-build/
CC=/usr/bin/gcc-4.2 CXX=/usr/bin/g++-4.2 ccmake ../BRAINSTools \
-DUSE_BRAINSConstellationDetector:BOOL=ON \
-DUSE_BRAINSABC:BOOL=ON
## NOTE: If you are using a version of Python different from the system default, CMake will ignore your
##       environment variables.  To ensure CMake points to the right version, you need to set the
##       following variables as well:
## -DPYTHON_LIBRARY:PATH=/path/to/python/lib \
## -DPYTHON_INCLUDE_DIR:PATH=/path/to/python/include
make -j24 -k ;
## NOTE: The fetching of data still has problems with parallel builds, so we need to restart it at least
#        once
make

## Testing ##
BRAINSTools_MAX_TEST_LEVEL adjusts how agressive the test suite is
so that long running tests or incomplete tests can easily be
silenced
1 - Run the absolute minimum very fast tests (These should always pass before any code commit)
3 - Run fast tests on continous builds (These need immediate attention if they begin to fail)
5 - Run moderate nightly tests (These need immediate attention if they begin to fail)
7 - Run long running extensive test that are a burden to normal development (perhaps test 1x per week)
8 - Run tests that fail due to incomplete test building, these are good ideas for test that we don't have time to make robust)
9 - Run silly tests that don't have much untility
set(BRAINSTools_MAX_TEST_LEVEL 3 CACHE STRING "Testing level for managing test burden")
