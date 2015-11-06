#!/bin/sh

export X3DNA=$HOME/programs/x3dna-v2.1/
x3dna_bin=$X3DNA/bin
find_pair=$x3dna_bin/find_pair
analyze=$x3dna_bin/analyze

$find_pair circ.pdb  bpfile_tmp.dat
$analyze bpfile_tmp.dat
