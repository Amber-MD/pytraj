#!/bin/sh


sh devtools/travis-ci/pyflakes_check.sh

PLATFORM=`python -c 'import sys; print(sys.platform)'`
echo "PLATFORM =" $PLATFORM

# set env
PYTRAJ_HOME=`pwd`
export CPPTRAJHOME=`pwd`"/cpptraj/"
cd tests/fake_amberhome/
export AMBERHOME=`pwd`

# go back to tests folder
cd ..

# print info
python run_simple_test.py

# run tests
if [ "$PLATFORM" = "linux" ]; then
    export LD_LIBRARY_PATH=`pwd`"/cpptraj/lib":$LD_LIBRARY_PATH
    nosetests --with-coverage --cover-package pytraj -vs test_*.py */test_*.py
    unset CPPTRAJHOME # for pytraj.tools coverage
    nosetests --with-coverage --cover-package pytraj -vs test_docs.py
else
    # osx (window too?), minimal tests
    nosetests -vs test_*.py */test_*.py
fi

# run examples
cd ../examples
python ./run_examples.py
cd mpi_examples

# only do mpi test for linux
if [ "$PLATFORM" = "linux" ]; then
    sh run_mpi_examples.sh 4 
fi

# go back to pytraj home
cd $PYTRAJ_HOME
