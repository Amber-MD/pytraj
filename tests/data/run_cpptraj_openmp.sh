#!/bin/sh

for nthreads in 1 2 4 8; do

export OMP_NUM_THREADS=$nthreads
echo $OMP_NUM_THREADS

cat >test.in <<EOF
parm Tc5b.top
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
trajin md1_prod.Tc5b.x
rms2d out test_openmp.Tc5b.n_threads_$nthreads.dat
EOF

$AMBERHOME/bin/cpptraj.OMP -i test.in >& log.openmp.$nthreads
done
