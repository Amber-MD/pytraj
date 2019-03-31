#!/bin/sh

git clone https://github.com/hainm/ambermini_test
mv ambermini_test ./tests/energies/


if [ "$PYPY" = "true" ]; then
    cwd=`pwd`
    fn=pypy-5.9-linux_x86_64-portable
    bz2fn=$fn.tar.bz2
    wget https://bitbucket.org/squeaky/portable-pypy/downloads/$bz2fn
    mkdir $HOME/pypy/
    mv $bz2fn $HOME/pypy
    cd $HOME/pypy
    tar -xf $bz2fn
    export PATH=$HOME/pypy/$fn/bin:$PATH 
    wget wget https://bootstrap.pypa.io/get-pip.py
    pypy get-pip.py
    pypy -m pip install cython==0.29 numpy
    pypy -m pip install pytest
    cd $cwd
else
    miniconda_version=4.3.21
    # download ambermini_test for energy calculation
    
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
    pip install -r pip-requirements.txt
    
    # create myenv
    conda create -y -n myenv python=$PYTHON_VERSION
    source activate myenv
    conda install -y --file conda-requirements.txt
    pip install -r pip-requirements.txt
fi
