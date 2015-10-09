#!/bin/sh

# PYTRAJHOME is the root folder of `pytraj`
export PYTRAJHOME=`pwd`
cd cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib

# turn off openmp. need to install pytraj with opnmp too. Too complicated.
make libcpptraj -j8 || bash ./configure -shared gnu || bash ./configure -amberlib -shared gnu || bash ./configure -nomathlib -shared gnu || make libcpptraj -j8 || exit 1
cd $PYTRAJHOME
