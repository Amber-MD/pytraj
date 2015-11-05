#!/bin/sh

# should run this in $PYTRAJHOME
# always install openmp as default
# User need to install libcpptraj manually if they want to disable openmp

install_type=$1
echo $install_type

# PYTRAJHOME is the root folder of `pytraj`
export PYTRAJHOME=`pwd`

if [[ $install_type == 'github' ]]; then
    git clone https://github.com/Amber-MD/cpptraj
    echo 'he he, installing from github'
fi

cd cpptraj/
export CPPTRAJHOME=`pwd`
cd $CPPTRAJHOME
mkdir lib

echo `pwd`

openmp=`python ../installs/check_openmp.py`
echo $openmp
# turn off openmp. need to install pytraj with openmp too. Too complicated.
bash ./configure -shared $openmp gnu || bash ./configure -amberlib -shared $openmp gnu ||
bash ./configure -nomathlib -shared $openmp gnu || exit 1

make libcpptraj -j8 || exit 1
cd $PYTRAJHOME

echo
echo
echo "make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
echo "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
echo "then 'python ./setup.py install'"
