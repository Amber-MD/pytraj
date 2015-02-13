#!/bin/sh

# PYTRAJHOME is the root folder of `pytraj`
export PYTRAJHOME=`pwd`
git clone http://github.com/mojyt/cpptraj
cd cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib
bash ./configure -nomathlib -shared gnu
make libcpptraj
cd $PYTRAJHOME

echo
echo
echo "make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
echo "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
echo "then 'python ./setup.py install'"
