#!/bin/bash

mkdir ./python_install_stuff
cd ./python_install_stuff

export INSTALL_DIR=/scratch/PREDICT/Experiments/NewExperiment/python-site-packages
export PYTHONPATH=${INSTALL_DIR}
export THIS_DIR=$(pwd)

## On 10.7.3 the gcc compiler is not reliable, use clang
#export CC=/usr/bin/clang
#export CXX=/usr/bin/clang++

if [ 1 -eq 1 ]; then
git clone http://github.com/numpy/numpy.git numpy
cd ${THIS_DIR}/numpy
python setup.py install --prefix=${INSTALL_DIR}
fi

if [ 0 -eq 1 ]; then
git clone http://github.com/scipy/scipy.git scipy
cd ${THIS_DIR}/scipy
python setup.py install --prefix=${INSTALL_DIR}
fi

for required_package in nibabel networkx nipype; do
  easy_install --install-dir=${INSTALL_DIR} ${required_package}
done

#git clone git://github.com/nipy/nipype.git  nipype
#cd ${THIS_DIR}/nipype
#python setup.py install --prefix=${INSTALL_DIR}


