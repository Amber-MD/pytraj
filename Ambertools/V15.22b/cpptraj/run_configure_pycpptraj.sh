#!/bin/sh

export AMBERHOME=/mnt/raidc/haichit/AMBER14_official.naga84.forPythonTest.2/
export LD_LIBRARY_PATH=$AMBERHOME/lib
./configure_pycpptraj -shared -nomathlib -amberlib gnu
