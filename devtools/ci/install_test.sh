#!/bin/sh

# for testing

if [ "$PYPY" = "true" ]; then
    echo "Nothing"
else
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
    if [ $osname = "darwin" ]; then
        brew install open-mpi
        pip install mpi4py
    else
        conda install mpi4py --y
    fi
fi
