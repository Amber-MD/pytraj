#!bin/sh

# create this file to hide output
# python setup.py install --amber-release


if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    python setup.py install --disable-openmp
else
    if [ "$TEST_SETUP" == 'true' ]; then
        echo "TEST_SETUP"
    else
        if [ "$USE_OPENMP" == 'false' ]; then 
            python setup.py install --disable-openmp
        else
            # pytraj will pick compiler based on COMPILER env
            # or using default (gnu for linux, clang for osx)
            python setup.py install
        fi
    fi
fi
