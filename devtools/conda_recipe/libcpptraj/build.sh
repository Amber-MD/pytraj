#!/bin/sh

bash ./configure -shared -nomathlib --with-netcdf=$PREFIX gnu
make libcpptraj

mkdir -p $PREFIX/lib/cpptraj/lib/
mkdir -p $PREFIX/lib/cpptraj/src/

cp lib/libcpptraj* $PREFIX/lib/cpptraj/lib/
cp lib/libcpptraj* $PREFIX/lib/

cp src/*.h $PREFIX/lib/cpptraj/src/
