#!/bin/sh

cython_version=0.26
# download ambermini_test for energy calculation
git clone https://github.com/hainm/ambermini_test
mv ambermini_test ./tests/energies/

osname=`python -c 'import sys; print(sys.platform)'`
if [ $osname = "darwin" ]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b

export PATH=$HOME/miniconda3/bin:$PATH
conda install --yes conda-build jinja2 anaconda-client pip numpy nomkl
pip install auditwheel
pip install cython==$cython_version

# create myenv
conda create -y -n myenv python=$PYTHON_VERSION
source activate myenv
conda update -y conda
conda install -y numpy h5py libnetcdf pyflakes
conda install -y hdf4 # temporary install to fix error "libmfhdf.so.0: cannot open shared object file"
pip install cython==$cython_version
pip install auditwheel
pip install mock
