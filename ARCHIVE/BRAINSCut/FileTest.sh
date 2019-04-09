#!/bin/bash
SuitePath=`pwd`
NO_ARGS=0

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    echo "Usage: `basename $0` options:"
    echo "-b    Binary Path"
    echo "-p    Build Path"
    echo "-d    Data Path"
    exit $E_OPTERROR        # Exit and explain usage, if no argument(s) given.
fi
while getopts "d:b:p:" Option
do
    case $Option in
  b ) BinaryPath=$OPTARG;;
        d ) DataPath=$OPTARG;;
        p ) BuildPath=$OPTARG;;
        * ) echo "Error: $0"
            echo "Unimplemented option chosen. Valid options are:"
            echo "-b    Binary Path"
            echo "-p    Build Path"
            echo "-d    Data Path"
            exit $E_OPTERROR;;
    esac
done


if [ ! -f ${BinaryPath}/VectorImageCompare  ]; then
    echo "FAILED to find ${BinaryPath}/VectorImageCompare"
    exit 1
fi

for x in Forward1.nii.gz Forward2.nii.gz Reverse1.nii.gz Reverse2.nii.gz
do
  #             VectorImageCompore
    if ${BinaryPath}/VectorImageCompare ${DataPath}/$x ${BuildPath}/Test$x
    then
        echo files match: ${DataPath}/$x ${BuildPath}/Test$x
    else
        echo files differ: ${DataPath}/$x ${BuildPath}/Test$x
        exit 1
    fi
done

for x in \
    Model.txtANNl_hippo000000010 \
    Model.txtANNr_hippo000000005 \
    Model.txtANNr_hippo000000020 \
    Model.txtANNl_hippo000000015 \
    Model.txtANNr_hippo000000010 \
    Model.txtANNl_hippo000000005 \
    Model.txtANNl_hippo000000020 \
    Model.txtANNr_hippo000000015
do
    if ${BinaryPath}/FloatListCompare ${DataPath}/$x ${BuildPath}/Test$x
    then
        echo files match: ${DataPath}/$x ${BuildPath}/Test$x
    else
        echo files differ: ${DataPath}/$x ${BuildPath}/Test$x
        exit 1
    fi
done
exit 0
