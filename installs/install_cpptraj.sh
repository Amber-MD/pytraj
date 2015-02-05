#!/bin/sh

# PYTRAJHOME is the root folder of `pytraj`

export PYTRAJHOME=`pwd`
#cd ./Ambertools/V15.22b/cpptraj/
cd ./Ambertools/dev/cpptraj
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib
#bash $PYTRAJHOME/installs/configure_pytraj -nomathlib -shared gnu
#bash $PYTRAJHOME/installs/cpptraj.V15.22b/configure_pytraj -nomathlib -shared -amberlib gnu
bash ./configure -nomathlib -shared -amberlib gnu
cd $CPPTRAJHOME/src
#make -f $PYTRAJHOME/installs/cpptraj.V15.22b/Makefile_libcpptraj libcpptraj
make libcpptraj
cd $PYTRAJHOME

echo
echo
echo "make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
echo "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
echo "then 'python ./setup.py install'"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CPPTRAJHOME/lib

# now try:
# python ./setup.py install

# if you got error: netcdf.h: No such file or directory
# you need to install netcdf lib, or specify the include dir, or usee "-amberlib"
