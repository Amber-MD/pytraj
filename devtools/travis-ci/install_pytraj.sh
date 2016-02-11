#!bin/sh

# create this file to hide output
# python setup.py install --amber-release

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    # git clone https://github.com/Amber-MD/cpptraj
    git clone https://github.com/swails/cpptraj
    cd cpptraj
    git checkout mactravis
    ./configure $BUILD_FLAGS clang
    make libcpptraj -j4
    cd ../
    python setup.py install --disable-openmp
else
    python setup.py install
fi
