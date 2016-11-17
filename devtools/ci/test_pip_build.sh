#!/bin/sh

./devtools/mkrelease
git clone https://github.com/amber-md/cpptraj
export CPPTRAJHOME=`pwd`/cpptraj/
conda install -y  \
    zlib \
    bzip2 \
    libnetcdf \
    openblas \
    libgfortran \
    gcc
pip install auditwheel
python scripts/install_libcpptraj.py -openmp
(cd dist && python ../scripts/build_wheel.py *gz)
