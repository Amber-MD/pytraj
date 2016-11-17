#!/bin/sh

./devtools/mkrelease
git clone https://github.com/amber-md/cpptraj
python scripts/install_libcpptraj.py -openmp
export CPPTRAJHOME=`pwd`/cpptraj/
(cd dist && python ../scripts/build_wheel.py *gz)
