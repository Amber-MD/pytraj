#!/bin/sh

if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-MacOSX-x86_64.sh -O miniconda.sh;
else
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
# install stable version
pip install conda

conda install --yes conda-build jinja2 anaconda-client pip

# create myenv
conda create -y -n myenv python=$PYTHON_VERSION numpy cython h5py mpi4py libnetcdf

source activate myenv
conda install --yes anaconda-client coverage pyflakes

if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    conda install netCDF4 -y
    conda update libnetcdf -y
    conda update netCDF4 -y
fi

# install other packages here
# conda install mdtraj -c omnia --yes
# pip install coveralls
# pip install coverage
# pip install nose
# pip install git+git://github.com/ParmEd/ParmEd
# pip install memory_profiler
# pip install psutil
# pip install cclib
