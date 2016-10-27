#!/bin/sh

PLATFORM=`uname -s | awk '{print $1}'`

if [[ $PLATFORM = "Darwin" ]] ; then
    bash configure -shared -macAccelerate --with-netcdf=/usr/local -noarpack clang
else
    bash ./configure --with-netcdf=$PREFIX \
                     --with-blas=$PREFIX \
                     --with-bzlib=$PREFIX \
                     --with-zlib=$PREFIX \
                     -openmp \
                     -shared \
                     -openblas -noarpack gnu
fi
make libcpptraj -j4

mkdir -p $PREFIX/include/cpptraj/

if [[ $PLATFORM = "Darwin" ]] ; then
    install_name_tool -id @rpath/libcpptraj.dylib lib/libcpptraj.dylib
    install_name_tool -change /usr/local/opt/netcdf/lib/libnetcdf.7.dylib @rpath/libnetcdf.7.dylib lib/libcpptraj.dylib
fi
cp lib/libcpptraj* $PREFIX/lib/
cp src/*.h $PREFIX/include/cpptraj/
