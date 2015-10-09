#!/bin/sh

# PYTRAJHOME is the root folder of `pytraj`
export PYTRAJHOME=`pwd`
git clone https://github.com/Amber-MD/cpptraj
cd cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib

# turn off openmp. need to install pytraj with openmp too. Too complicated.
bash ./configure -shared gnu || bash ./configure -amberlib -shared gnu || bash ./configure -nomathlib -shared gnu || exit 1
make libcpptraj -j8 || exit 1
cd $PYTRAJHOME

echo
echo
echo "make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
echo "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
echo "then 'python ./setup.py install'"
