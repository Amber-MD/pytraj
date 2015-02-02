#!/bin/sh

#wrapp vim for openning two files

pxd=$1
cpptrajsrc="/mnt/raidc/haichit/AMBER14_official.naga84.forPythonTest.2/AmberTools/src/cpptraj/src/"

if [ ! -f ${pxd}.pxd  ]; then
    cat PXD_template.dat > ${pxd}.pxd
fi

vim -O ${pxd}.pxd  $cpptrajsrc${pxd}.h
