#!bin/sh

# create this file to hide output
# python setup.py install --amber-release

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    python setup.py install --disable-openmp
else
    if [ "$COMPILER" == "clang" ]; then
        CC=clang CXX=clang++ python setup.py install --disable-openmp
    else
        python setup.py install
    fi
fi
