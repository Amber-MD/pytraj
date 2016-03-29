#!/bin/sh

#wrap vim for opening two files

#rootname=$1
#cpptrajsrc=$AMBERHOME"/AmberTools/src/cpptraj/src/"
cpptrajsrc=$CPPTRAJHOME/src/

#if [ ! -f ${rootname}.pyx ]; then
#  cat PYX_template.dat > ${rootname}.pyx
#fi

#vim -O ${rootname}.pyx  $cpptrajsrc${rootname}.cpp
pyx=$1
cpptrajfile=$cpptrajsrc$2
vim -O $pyx $cpptrajfile
