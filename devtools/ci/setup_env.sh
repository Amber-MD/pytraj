#!/bin/sh

cython_version=0.26
miniconda_version=4.3.21
# download ambermini_test for energy calculation
git clone https://github.com/hainm/ambermini_test
mv ambermini_test ./tests/energies/

osname=`python -c 'import sys; print(sys.platform)'`
if [ $osname = "darwin" ]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-${miniconda_version}-MacOSX-x86_64.sh -O miniconda.sh;
else
    wget https://repo.continuum.io/miniconda/Miniconda3-${miniconda_version}-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b

export PATH=$HOME/miniconda3/bin:$PATH
conda install --yes conda-build=3.0.19 jinja2 anaconda-client pip numpy=1.13.1 nomkl
pip install auditwheel==1.7.0
pip install cython==$cython_version

# create myenv
conda create -y -n myenv python=$PYTHON_VERSION
source activate myenv
# conda install -y numpy=1.13.1 nomkl h5py=2.7.0 libnetcdf=4.4.1  pyflakes
conda install -y conda=4.3.25
conda install -y libnetcdf=4.4.1
pip install pyflakes==1.6.0
pip install numpy==1.13.1
pip install h5py==2.7.0
conda install -y hdf4=4.2.12 # temporary install to fix error "libmfhdf.so.0: cannot open shared object file"
pip install cython==$cython_version
pip install auditwheel==1.7.0
pip install mock==2.0.0
