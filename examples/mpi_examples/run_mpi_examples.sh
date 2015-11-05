#!/bin/sh

cores=$1
mpirun -n $cores python mpi_calc_molsurf_0_single_traj.py || exit 1
mpirun -n $cores python mpi_load_batch.py || exit 1
mpirun -n $cores python mpi_search_hbonds_single_traj.py || exit 1
mpirun -n $cores python mpi_pysander.py || exit 1
echo "mpi ok"
