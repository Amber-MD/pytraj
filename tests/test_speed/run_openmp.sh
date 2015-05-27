#!/bin/sh

my_script=$1

for n_cpus in 1 2 4 6 8; do
    echo "n_cpus = " $n_cpus
    export OMP_NUM_THREADS=$n_cpus
    python $my_script
done
