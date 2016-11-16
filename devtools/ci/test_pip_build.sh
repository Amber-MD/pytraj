#!/bin/sh

./devtools/mkrelease
git clone https://github.com/amber-md/cpptraj
export CPPTRAJHOME=`pwd`/cpptraj/
python scripts/install_libcpptraj.py -openmp
(cd dist && python ../scripts/build_wheel.py *gz)
