#!/bin/sh

# download ambermini_test for energy calculation
git clone https://github.com/hainm/ambermini_test
mv ambermini_test ./tests/energies/

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
ls $HOME/miniconda/bin/
conda install --yes conda-build jinja2 anaconda-client pip
conda install -yes numpy cython=0.26 h5py libnetcdf
source devtools/ci/install_test.sh
