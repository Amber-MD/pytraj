#!/bin/sh

# to enable openmp in pytraj, just need to compile `libcpptraj` with openmp
# (that's it)
# cd $AMBERHOME
# ./configure -openmp gnu 
# cd $AMBERHOME/AmberTools/src
# make pytraj

for nthreads in 1 2 4 8; do
    export OMP_NUM_THREADS=$nthreads
    echo "n_threads = " $OMP_NUM_THREADS
    python ./test_openmp_0.py
done
