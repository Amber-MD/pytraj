#!/bin/sh

wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
# install stable version
pip install conda

conda install --yes conda-build jinja2 anaconda-client pip

# create myenv
conda create -y -n myenv python=$PYTHON_VERSION numpy cython h5py mpi4py libnetcdf

source activate myenv
conda install --yes anaconda-client coverage pyflakes
conda install mdtraj -c omnia --yes

# install other packages here
pip install coveralls
pip install coverage
pip install nose
pip install git+git://github.com/ParmEd/ParmEd
pip install memory_profiler
pip install psutil
pip install cclib
