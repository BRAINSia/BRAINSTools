BRAINSTools build configuration scripts
=======================================


Helpful build starting points for common deployment and testing environments.

```
mkdir -p src
cd src
git clone git@github.com:BRAINSia/BRAINSTools.git
mkdir -p BRAINSTools-llvm10-RelWithDebInfo
cmake \
    -G Ninja \
    -S BRAINSTools \
    -B BRAINSTools-RelWithDebInfo \
    -C BRAINSTools/Utilities/build_configs/Darwin/BRAINSTools-llvm10-RelWithDebInfo.cmake

cd BRAINSTools-llvm10-RelWithDebInfo
ninja

```
