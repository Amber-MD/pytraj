#!/bin/sh

./devtools/mkrelease
wget https://github.com/amber-md/cpptraj/archive/master.zip
unzip master.zip
export CPPTRAJHOME=`pwd`/cpptraj/
python scripts/install_libcpptraj.py -openmp
(cd dist && python ../scripts/build_wheel.py *gz)
