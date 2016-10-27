#!/bin/sh

# download ambermini_test for energy calculation
git clone https://github.com/hainm/ambermini_test
mv ambermini_test ./tests/energies/

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
ls $HOME/miniconda/bin/
conda install --yes conda-build jinja2 anaconda-client pip
conda install -yes numpy cython h5py libnetcdf
conda install --yes anaconda-client coverage pyflakes jupyter notebook

# install other packages here
conda install parmed -c ambermd --yes
conda install pysander -c hainm --yes
conda install cclib -c omnia --yes
conda install nglview -c bioconda --yes
pip install coveralls
pip install coverage
pip install pytest-cov
pip install nose
pip install memory_profiler
pip install psutil
# pip install cclib # install error with current v1.5 (08/24/2016)
pip install tqdm

if [ "$TRAVIS_OS_NAME" = "linux" ]; then
    # only test mpi on linux
    conda install mpi4py --yes
else
    # osx
    brew install open-mpi
    pip install mpi4py
fi
