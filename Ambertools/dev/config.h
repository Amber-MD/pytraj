# config.h for cpptraj
# configured using: ./configure -nomathlib -shared -amberlib gnu

CPPTRAJHOME=/mnt/raidc2/haichit/Study/Cython/pytraj_git_fork/pytraj/Ambertools/dev
CPPTRAJBIN=/mnt/raidc2/haichit/Study/Cython/pytraj_git_fork/pytraj/Ambertools/dev/bin
CPPTRAJLIB=/mnt/raidc2/haichit/Study/Cython/pytraj_git_fork/pytraj/Ambertools/dev/lib

DBGFLAGS=
CC=gcc
CXX=g++
FC=gfortran
CFLAGS= -O3 -Wall  -DNO_MATHLIB -DHASBZ2 -DHASGZ -DBINTRAJ -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/mnt/raidc/haichit/AMBER14_official.naga84.forPythonTest.2//include -fPIC $(DBGFLAGS)
CXXFLAGS= -O3 -Wall  -DNO_MATHLIB -DHASBZ2 -DHASGZ -DBINTRAJ -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/mnt/raidc/haichit/AMBER14_official.naga84.forPythonTest.2//include -fPIC $(DBGFLAGS)
FFLAGS=-ffree-form  -O3  -DNO_MATHLIB -DHASBZ2 -DHASGZ -DBINTRAJ -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/mnt/raidc/haichit/AMBER14_official.naga84.forPythonTest.2//include -fPIC $(DBGFLAGS)

READLINE=readline/libreadline.a
READLINE_HOME=readline

LDFLAGS=    -L/mnt/raidc/haichit/AMBER14_official.naga84.forPythonTest.2//lib -lnetcdf  -lbz2 -lz  -lgfortran -w
SFX=
