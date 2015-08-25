The BRAINSTools is a harness to assist in building the many of the BRAINSTools under development.

Developers should run the `./Utilities/SetupForDevelopment.sh` script to get started.

If developing on Mac OS X, make sure the xcode command line tools are installed with the command:
`xcode-select --install`

For more information on the individual BRAINSTools, please see the following link:
https://github.com/BRAINSia/BRAINSTools/wiki

## Building
### Mac OSX
Building BRAINSTools on mac is the same as any standard out of
source cmake build. Make a folder for the BRAINSTools project.
Clone the BRAINSTools repository, and run the set-up script. Make
a build directory contained in the same folder as the git repository.
From the build directory, run cmake (or ccmake to manually configure
the build settings) with the BRAINSTools repo as an argument, then
run make. For Example:

```sh
mkdir brains
cd brains
git clone http://github.com/BRAINSIa/BRAINSTools.git

cd BRAINSTools
./Utilities/SetupForDevelopment.sh
cd ..

mkdir build
cd build
cmake ../BRAINSTools

## NOTE: To configure options in the cmake GUI, use:
##     ccmake ../BRAINSTools
## instead
make -j${NUMOFTHREADS} -k

## NOTE: To find the number of threads from the OSX terminal, use:
##    sysctl -n hw.ncpu
```

### Linux RedHat 6
Example session for a clean build:

```sh
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

make -j${NUMOFTHREADS} -k

## NOTE: The fetching of data still has problems with parallel builds, so we need to restart it at least
#        once
make
```

### Linux Debian (Ubuntu)
Building BRAINSTools on a fresh install has the additional dependency
of building CMake on your system.  You cannot use the version from
`apt-get` as that does some unnatural things with Python resources to
be backwards compatible (see http://public.kitware.com/Bug/view.php?id=14156).

1) Install the necessary dependencies:
```sh
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install git python2.7
python2.7-dev g++ freeglut3-dev
```

2) Get the CMake binaries:
```sh
# You can find the URL of the latest CMake binary by
# * Go to www.cmake.org/download in a browser
# * Hover over the desired 'Binary distribution' for your operating
#     system
# * Right-click and select 'Copy Link Address' for the file ending in
#     "tar.gz"
URL=http://www.cmake.org/files/v3.2/cmake-3.2.1-Linux-x86_64.tar.gz
wget ${URL}
# Decompress the file
tar -xzvf cmake-3.2.1-Linux-x86_64.tar.gz
# Add the binary to your PATH environment variable
export PATH=${PWD}/cmake-3.2.1-Linux-x86_64/bin:${PATH}
```

2) Clone the repository and build:
```sh
git clone https://github.com/BRAINSia/BRAINSTools.git
mkdir build
cd build
CC=/usr/bin/gcc-4.8 \
CXX=/usr/bin/g++-4.8 \
cmake ../BRAINSTools \
make -j${NUMOFTHREADS} -k
```
:warning: You can find the number of threads on your system in Ubuntu with `lscpu`

## Testing
`BRAINSTools_MAX_TEST_LEVEL` adjusts how agressive the test suite is
so that long running tests or incomplete tests can easily be
silenced

```
1 - Run the absolute minimum very fast tests (These should always pass before any code commit)
3 - Run fast tests on continous builds (These need immediate attention if they begin to fail)
5 - Run moderate nightly tests (These need immediate attention if they begin to fail)
7 - Run long running extensive test that are a burden to normal development (perhaps test 1x per week)
8 - Run tests that fail due to incomplete test building, these are good ideas for test that we don't have time to make robust)
9 - Run silly tests that don't have much untility
set(BRAINSTools_MAX_TEST_LEVEL 3 CACHE STRING "Testing level for managing test burden")
```
