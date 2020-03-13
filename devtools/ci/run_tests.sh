#!/bin/sh


export OMP_NUM_THREADS=1

if [ "$TEST_SETUP" = 'true' ]; then
    echo "Test setup command line"
    py.test -vs devtools/ci/test_setup_command.py
    # run this on circleci
    # sh devtools/ci/test_pip_build.sh
else
    if [ "$PYPY" = "true" ]; then
        cd tests && pypy -c "import pytraj; pytraj.run_tests()"
    else
        sh devtools/ci/pyflakes_check.sh || exit 1
        
        isOSX=`python -c 'import sys; print(sys.platform.startswith("darwin"))'`
        echo "isOSX =" $isOSX
        
        # set env
        PYTRAJ_HOME=`pwd`
        export CPPTRAJHOME=`pwd`"/cpptraj/"
        cd tests/fake_amberhome/
        export AMBERHOME=`pwd`
        
        cd $PYTRAJ_HOME
        if [ "$isOSX" = "True" ]; then
            echo "Minimal tests for OSX"
            (cd tests && py.test -v test_analysis/)
        else
            python run_tests.py -c || exit 1
        fi
        
        # run examples
        cd ./examples
        python ./run_examples.py || exit 1
        cd mpi_examples
        
        # only do mpi test for linux
        if [ "$isOSX" != "True" ]; then
            sh run_mpi_examples.sh 1 || exit 1
        fi
        
        # go back to pytraj home
        cd $PYTRAJ_HOME

        # test pip install on osx
        if [ "$isOSX" = "True" ]; then
            echo "Testing wheel building"
            # back to root env
            # (seems conda-build like root env)
            source deactivate
            ./devtools/mkrelease
            (cd dist && python ../scripts/build_wheel.py ./pytraj*gz)
        fi
    fi
fi
