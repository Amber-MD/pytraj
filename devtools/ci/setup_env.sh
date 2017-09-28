#!/bin/sh

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
conda install -y --file conda-requirements.txt
pip install -r pip-requirements.txt
