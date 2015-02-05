#!/bin/sh

# PYTRAJHOME is the root folder of `pytraj`
export PYTRAJHOME=`pwd`
cd ./Ambertools/dev/cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib
bash ./configure -nomathlib -shared gnu
cd $CPPTRAJHOME/src
make libcpptraj
cd $PYTRAJHOME

echo
echo
echo "make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
echo "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
echo "then 'python ./setup.py install'"
