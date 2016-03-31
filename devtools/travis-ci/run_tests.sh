#!/bin/sh

python -c 'import pytraj as pt'
python ./run_simple_test.py

if [ "$TRAVIS_OS_NAME" = "linux" ]; then
    nosetests --with-coverage --cover-package pytraj -vs test_*.py */test_*.py
    unset CPPTRAJHOME # for pytraj.tools coverage
    nosetests --with-coverage --cover-package pytraj -vs test_docs.py
    cd ../examples
    python ./run_examples.py
    cd mpi_examples
    sh run_mpi_examples.sh 4 
else
    # osx (window too?), minimal tests
    echo "osx test"
    nosetests -vs .
