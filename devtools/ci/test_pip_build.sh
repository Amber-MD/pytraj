#!/bin/sh

./devtools/mkrelease
export CPPTRAJHOME=`pwd`/cpptraj/
python scripts/install_libcpptraj.py -openmp
(cd dist && python ../scripts/build_wheel.py *gz)
