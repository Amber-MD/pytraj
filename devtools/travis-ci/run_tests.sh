#!/bin/sh


sh devtools/travis-ci/pyflakes_check.sh || exit 1

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
python run_simple_test.py || exit 1

# run tests
if [ "$PLATFORM" = "linux" ]; then
    nosetests --with-coverage --cover-package pytraj -vs test_*.py */test_*.py || exit 1
    unset CPPTRAJHOME # for pytraj.tools coverage
    nosetests --with-coverage --cover-package pytraj -vs test_docs.py || exit 1
else
    # osx (window too?), minimal tests
    nosetests -vs test_*.py */test_*.py || exit 1
fi

# run examples
cd ../examples
python ./run_examples.py || exit 1
cd mpi_examples

# only do mpi test for linux
if [ "$PLATFORM" = "linux" ]; then
    sh run_mpi_examples.sh 4 || exit 1
fi

# go back to pytraj home
cd $PYTRAJ_HOME
