#!/bin/sh

# for testing
conda install -y anaconda-client jupyter notebook
conda install -y parmed -c ambermd
conda install -y pysander -c hainm
conda install -y cclib -c omnia
conda install -y nglview -c bioconda
pip install coveralls
pip install coverage
pip install pytest-cov
pip install nose
pip install memory_profiler
pip install psutil
pip install tqdm

osname=`python -c 'import sys; print(sys.platform)'`
if [ $osname = "linux" ]; then
    # only test mpi on linux
    # conda install mpi4py --y
    pip install mpi4py
else
    # osx
    brew install open-mpi
    pip install mpi4py
fi
