#!bin/sh

# create this file to hide output
# python setup.py install --amber-release

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    python setup.py install --disable-openmp
else
    if [[ "$TEST_SETUP" == 'true' ]]; then
        echo "TEST_SETUP"
    else
        if [ "$COMPILER" == "clang" ]; then
            export CC=clang
            export CXX=clang++
            python setup.py install --disable-openmp
        else
            python setup.py install
        fi
    fi
fi
