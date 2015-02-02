#!/bin/sh

# PYCPPTRAJ_HOME is the root folder of `pycpptraj`

export PYCPPTRAJ_HOME=`pwd`
cd ./Ambertools/V15.22b/cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib
#bash $PYCPPTRAJ_HOME/installs/configure_pycpptraj -nomathlib -shared gnu
bash $PYCPPTRAJ_HOME/installs/cpptraj.V15.22b/configure_pycpptraj -nomathlib -shared -amberlib gnu
cd $CPPTRAJHOME/src
make -f $PYCPPTRAJ_HOME/installs/cpptraj.V15.22b/Makefile_libcpptraj libcpptraj
cd $PYCPPTRAJ_HOME

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
