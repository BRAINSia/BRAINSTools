The BRAINSTools project is a harness to assist in building the many BRAINSTools under development.

## Building
Example session for a clean build on Linux using gcc-4.2:
```bash
DIRECTORY=~/src/mytester
mkdir ${DIRECTORY}
cd ${DIRECTORY}
git clone git://github.com/BRAINSia/BRAINSTools.git
```
__For developers:__

```bash
cd ${DIRECTORY}/BRAINSTools/
bash ./Utilities/SetupForDevelopment.sh
```

Finally:
```bash
mkdir -p ${DIRECTORY}/BRAINSTools-build
cd ${DIRECTORY}/BRAINSTools-build/
CC=/usr/bin/gcc-4.2 \
CXX=/usr/bin/g++-4.2 \
ccmake ../BRAINSTools \
-DUSE_BRAINSFit:BOOL=ON \
-DUSE_BRAINSConstellationDetector:BOOL=ON \
-DUSE_BRAINSABC:BOOL=ON
make -j24 -k
## NOTE:  The fetching of data still has problems 
## with parallel builds, so you may need to restart it a few times...
make
make
make
```

# Development
Developers should run the "./Utilities/SetupForDevelopment.sh" script to get started.

## Testing
`BRAINSTools_MAX_TEST_LEVEL` adjusts how agressive the test suite is so that long running tests or incomplete tests can easily be silenced

```cmake
set(BRAINSTools_MAX_TEST_LEVEL 3 
      CACHE STRING "Testing level for managing test burden")
```

__1__ - Run the absolute minimum (very fast tests) 
  * These should always pass before any code commit!

__3__ - Run fast tests on continous builds.
* These need immediate attention if they begin to fail!

__5__ - Run moderate nightly tests.
  * These need immediate attention if they begin to fail!

__7__ - Run long running extensive test that are a burden to normal development.
  * Testing done 1x per week.

__8__ - Run tests that fail due to incomplete test building. 
  * These are good ideas for test that we don't have time to make robust.

__9__ - Run tests that don't have much utility currently.

***

###### Example
setting a test's max level in TestSuite/CMakeLists.txt
```cmake
if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 8)
  ExternalData_add_test(FindCenterOfBrainFetchData
    NAME itkResampleInPlaceImageFilterTest
    COMMAND $<TARGET_FILE:itkResampleInPlaceImageFilterTest>
      itkResampleInPlaceImageFilterTest input1 transform1 checkresult
  )
```