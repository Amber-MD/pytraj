#!/bin/sh


export OMP_NUM_THREADS=1

if [ "$TEST_SETUP" == 'true' ]; then
    echo "Test setup command line"
    py.test -vs devtools/ci/test_setup_command.py
else
    sh devtools/ci/pyflakes_check.sh || exit 1
    
    PLATFORM=`python -c 'import sys; print(sys.platform)'`
    echo "PLATFORM =" $PLATFORM
    
    # set env
    PYTRAJ_HOME=`pwd`
    export CPPTRAJHOME=`pwd`"/cpptraj/"
    cd tests/energies/fake_amberhome/
    export AMBERHOME=`pwd`
    
    cd $PYTRAJ_HOME
    python run_tests.py -c || exit 1
    
    # run examples
    cd ./examples
    python ./run_examples.py || exit 1
    cd mpi_examples
    
    # only do mpi test for linux
    if [ "$PLATFORM" = "linux" ]; then
        sh run_mpi_examples.sh 4 || exit 1
    fi
    
    # go back to pytraj home
    cd $PYTRAJ_HOME
fi
