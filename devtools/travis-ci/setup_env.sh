#!/bin/sh

if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-MacOSX-x86_64.sh -O miniconda.sh;
else
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
conda install --yes conda-build jinja2 anaconda-client pip

# create myenv
conda create -y -n myenv python=$PYTHON_VERSION numpy cython h5py libnetcdf

source activate myenv
conda install --yes anaconda-client coverage pyflakes jupyter notebook

# install other packages here
conda install parmed -c ambermd --yes
conda install pysander -c hainm --yes
pip install git+https://github.com/arose/nglview
pip install coveralls
pip install coverage
pip install nose
pip install memory_profiler
pip install psutil
pip install cclib

if [ "$TRAVIS_OS_NAME" = "linux" ]; then
    # only test mpi on linux
    conda install mpi4py --yes
else
    # osx
    brew install open-mpi
    pip install mpi4py
fi
