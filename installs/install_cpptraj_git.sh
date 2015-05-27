#!/bin/sh

# PYTRAJHOME is the root folder of `pytraj`
export PYTRAJHOME=`pwd`
git clone https://github.com/mojyt/cpptraj
cd cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib
# TODO: smarter picking flags
#bash ./configure -shared -openmp gnu || bash ./configure -nomathlib -shared -openmp gnu \
#|| bash ./configure -nomathlib -shared -nonetcdf -openmp gnu || bash ./configure -shared gnu \
#|| bash ./configure -nomathlib -shared gnu || bash ./configure -nomathlib -shared -nonetcdf gnu

# turn off openmp. need to install pytraj with opnmp too. Too complicated.
bash ./configure -shared gnu \
|| bash ./configure -nomathlib -shared gnu || bash ./configure -nomathlib -shared -nonetcdf gnu
make libcpptraj -j8
cd $PYTRAJHOME

echo
echo
echo "make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
echo "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
echo "then 'python ./setup.py install'"
