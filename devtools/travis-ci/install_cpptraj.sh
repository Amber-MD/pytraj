#!/bin/sh

if [[ $CPPTRAJ_ANACONDA == "NO" ]]; then
    # build libcpptraj on travis
    conda build devtools/conda-recipe/libcpptraj
    conda install $HOME/miniconda/conda-bld/linux-64/libcpptraj-dev-* --yes
else
    # install from anaconda
    # this is for testingg lapack, arpack, blas
    conda install -c ambermd libcpptraj-dev --force --yes
fi
