#!/bin/sh

bash ./configure -shared -nomathlib --with-netcdf=$PREFIX -openmp gnu
make libcpptraj

mkdir -p $PREFIX/include/cpptraj/

cp lib/libcpptraj* $PREFIX/lib/

cp src/*.h $PREFIX/include/cpptraj/
