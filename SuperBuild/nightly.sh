#!/bin/bash
#
# Set the library paths such that all locally built shared 
# libraries are found and used
# ahead of system libs
package=BRAINSStandAlone

#
# when run by cron, the path variable is only /bin:/usr/bin
export PATH="/opt/cmake/bin:/usr/local/bin:/usr/sbin:$PATH"

#
# make the testing directory, based on current user name
#user=`who -m | sed -e 's/ .*$//'`
user=${LOGNAME}

ThisComputer=`hostname`


#
# the default is to use /brainsdev/kent -- which is
# appropriate on the b2dev VMs.
if [ $# = 0 ] ; then
    startdir=/brainsdev/kent/Testing
else
    startdir=$1
    shift
fi

#
# needed for ssl authentication for git
export GIT_SSL_NO_VERIFY=true

CXXFLAGS="${CXXFLAGS:-}"
CFLAGS="${CFLAGS:-}"
LDFLAGS="${LDFLAGS:-}"
CCOverride=""
CXXOverride=""

while [ $# -gt 0 ]
do
    echo Current Arg $1
    case "$1" in
        coverage)
            coverage=1 ;;
        CC=*)
          echo Override C Compiler ; CCOverride=`echo $1 | sed -e 's/^CC=//'` ;;
        CXX=*)
          echo Override CXX Compiler ; CXXOverride=`echo $1 | sed -e 's/^CXX=//'` ;;
    esac
    shift
done

OsName=$(uname)
if [ "${OsName}" = "Darwin" ] ; then
    Compiler=clang-`clang -v 2>&1 | head -1 | awk '{print $4}'`
    Compiler=${Compiler}-`clang -v 2>&1 | tail -2 | head -1 | awk '{print $2}'`
    export CC=`which clang`
    export CXX=`which clang++`
else
    which gcc > /dev/null 2>&1
    if [ $? == 0 ] ; then
        Compiler=gcc-`gcc -dumpversion`-`gcc -dumpmachine`
    else
        Compiler=unknown
    fi
fi

if [ X$CCOverride != X ] ; then
    export CC=$CCOverride;
    Compiler=`${CCOverride} -dumpversion`-`gcc -dumpmachine`
fi
if [ X$CXXOverride != X ] ; then
    export CXX=$CXXOverride;
    Compiler=`${CXXOverride} -dumpversion`-`gcc -dumpmachine`
fi

echo "Compiler=${Compiler} CC=${CC} CXX=${CXX} coverage=${coverage}"

OS=$(uname -s)
NPROCS=1

if [ "${OS}" = "Linux" ] ; then
#    NPROCS=$(grep -c ^processor /proc/cpuinfo)
    export CFLAGS="${CFLAGS} -fpic"
    export CXXFLAGS="${CXXFLAGS} -fpic"
#else
#    NPROCS=$(system_profiler | awk '/Number Of Cores/{print $5}{next;}')
fi

# create the testing directory if necessary
mkdir -p ${startdir}
if [ ! -d ${startdir} ] ; then
    echo ${startdir} cannot be created, exiting
    exit 1
fi

cd ${startdir}

echo checking out test data in `pwd`

mkdir -p ${startdir}/${ThisComputer}

cd ${startdir}/${ThisComputer}

top=`pwd`
echo WORKING IN $top

# check out BRAINSStandAlone in a directory unique to each host -- this is unfortunately necessary
# because svn can't update a directory  checked out by a newer version of svn, so
# every host has their own copy of BRAINS3 so that it's compatible with the local svn version.
if [ -d BRAINSStandAlone ] ; then
    cd BRAINSStandAlone
    git pull
else
    git clone git@github.com:BRAINSia/BRAINSStandAlone.git
fi
if [ $? != 0 ]
then
    echo BRAINSStandAlone checkout failed, continuing with old version
fi



for BUILD_TYPE in Debug Release
do
    B3Build=${top}/${BUILD_TYPE}-${Compiler}
    if [ "$BUILD_TYPE" = "Debug" -a "$coverage" = "1" ] ; then
	CXXFLAGS="${CXXFLAGS} -g -O0 -Wall -W -Wshadow -Wunused-variable \
	    -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers \
	    -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage"
	CFLAGS="${CFLAGS} -g -O0 -Wall -W -fprofile-arcs -ftest-coverage"
	LDFLAGS="${LDFLAGS} -fprofile-arcs -ftest-coverage"
    fi
    mkdir -p ${B3Build}
    cd ${B3Build}
    rm -f CMakeCache.txt
    # force reconfigure.
    find . -name '*-configure' | xargs rm -f
    #
    # the Build type
    cmake -DSITE:STRING=${ThisComputer} \
        -G "Unix Makefiles" \
        -DCMAKE_CXX_COMPILER:STRING="${CXX}" \
        -DCMAKE_C_COMPILER:STRING="${CC}" \
	-DCMAKE_C_FLAGS:STRING="${CFLAGS}" \
	-DCMAKE_CXX_FLAGS:STRING="${CXXFLAGS}" \
	-DCMAKE_EXE_LINKER_FLAGS:STRING="${LDFLAGS}" \
	-DCMAKE_MODULE_LINKER_FLAGS:STRING="${LDFLAGS}" \
	-DCMAKE_SHARED_LINKER_FLAGS:STRING="${LDFLAGS}" \
	-DBUILDNAME:STRING="${OsName}-${Compiler}-${BUILD_TYPE}" \
        -DBUILD_SHARED_LIBS:BOOL=Off \
	-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
        ${top}/BRAINSStandAlone
    echo "Building in `pwd`"
    scriptname=`basename $0`
    make -j ${NPROCS}
    cd BRAINSTools-build
    make clean
    if [ $scriptname = "nightly.sh" ] ; then
	ctest -j ${NPROCS} -D Nightly
    else
	ctest -j ${NPROCS} -D Experimental
    fi
    cd ..
done

cd ${top}
