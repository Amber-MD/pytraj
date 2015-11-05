#!/bin/sh

cores=$1
mpirun -n $cores python mpi_molsurf.py || exit 1
mpirun -n $cores python mpi_load_batch.py || exit 1
mpirun -n $cores python mpi_pysander.py || exit 1
echo "mpi ok"
