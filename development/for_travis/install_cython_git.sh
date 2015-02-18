#!/bin/sh

set -ex
export PYCPPTRAJ_HOME=`pwd`
git clone https://github.com/cython/cython 
cd cython*
python ./setup.py install
cd $PYCPPTRAJ_HOME
