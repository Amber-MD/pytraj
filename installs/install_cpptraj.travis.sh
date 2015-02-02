#!/bin/sh

# PYCPPTRAJ_HOME is the root folder of `pycpptraj`

export PYCPPTRAJ_HOME=`pwd`
cd ./Ambertools/v14/cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib
#bash $PYCPPTRAJ_HOME/installs/Ambertolls14/configure_pycpptraj -nomathlib -shared -amberlib gnu
bash $PYCPPTRAJ_HOME/installs/Ambertolls14/configure_pycpptraj -nomathlib -shared gnu
cd $CPPTRAJHOME/src
make -f $PYCPPTRAJ_HOME/installs/Ambertolls14/Makefile_libcpptraj libcpptraj
cd $PYCPPTRAJ_HOME

echo
echo
echo "make sure to export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$CPPTRAJHOME/lib"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CPPTRAJHOME/lib

# now try:
# python ./setup.py install

# if you got error: netcdf.h: No such file or directory
# you need to install netcdf lib, or specify the include dir, or usee "-amberlib"
