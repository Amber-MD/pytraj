#!/bin/sh

# for testing
conda install -y anaconda-client jupyter notebook
pip install cclib
pip install nglview
pip install coveralls
pip install coverage
pip install pytest-cov
pip install nose
pip install tqdm

source devtools/ci/install_ambertools.sh

osname=`python -c 'import sys; print(sys.platform)'`
if [ $osname = "linux" ]; then
    # only test mpi on linux
    conda install mpi4py --y
    # pip install mpi4py # compiling error (09/27/2017)
else
    # osx
    brew install open-mpi
    pip install mpi4py
fi
