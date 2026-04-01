# BRAINSTools Build System

## Standard SuperBuild

BRAINSTools requires CMake 3.20.6+ and a C++17 compiler (clang or gcc).
Builds must be out-of-source.

```bash
# Clone and set up
git clone https://github.com/BRAINSia/BRAINSTools.git
cd BRAINSTools
./Utilities/SetupForDevelopment.sh   # install git hooks; required once

# Configure (Phase I superbuild)
mkdir ../BRAINSTools-build
cd ../BRAINSTools-build
cmake ../BRAINSTools \
  -DCMAKE_BUILD_TYPE=Release \
  -DBRAINSTools_SUPERBUILD=ON

# Build (Phase I fetches deps, Phase II builds tools)
make -j$(nproc)
```

## Key CMake Options

```cmake
BRAINSTools_SUPERBUILD=ON       # Two-phase superbuild (default ON)
BRAINSTools_PKGBUILD=ON         # Inner product build (set automatically)

# System library overrides
USE_SYSTEM_ITK=OFF              # Use system ITK instead of building from source
USE_SYSTEM_VTK=OFF
USE_SYSTEM_SlicerExecutionModel=OFF
USE_SYSTEM_zlib=OFF

# Feature flags (set in Common.cmake)
USE_BRAINSFit=ON
USE_BRAINSABC=ON
USE_BRAINSConstellationDetector=ON
# ... see Common.cmake for full list

# Coverage build
BUILD_COVERAGE=ON               # Requires -DCMAKE_BUILD_TYPE=Debug
```

## Overriding a Dependency Source

To test against a local ITK fork without editing source:
```bash
cmake . \
  -DBRAINSTools_ITKv5_GIT_REPOSITORY=https://github.com/myuser/ITK.git \
  -DBRAINSTools_ITKv5_GIT_TAG=my-fix-branch
```
`ExternalProject_SetIfNotDefined` in each `SuperBuild/External_*.cmake` file
makes this possible. The variable name pattern is always
`BRAINSTools_<ProjName>_GIT_REPOSITORY` and `BRAINSTools_<ProjName>_GIT_TAG`.

## Building Only the Inner Product

After a full superbuild, iterate on BRAINSTools code without rebuilding deps:

```bash
cd BRAINSTools-build/BRAINSTools-Release-<version>-build
make -j$(nproc)
```

## Install Layout

```
BRAINSTools-build/BRAINSTools-Release-5.8.0/   ← CMAKE_INSTALL_PREFIX
├── bin/          ← compiled tool executables
├── lib/          ← shared libraries
│   └── cmake/    ← CMake config files
└── include/      ← public headers
```

## Build Config Files

Example platform-specific cache-init files are in
`Utilities/build_configs/Darwin/`. Pass with:
```bash
cmake -C ../BRAINSTools/Utilities/build_configs/Darwin/BRAINSTools-llvm10-Release.cmake ../BRAINSTools
```
