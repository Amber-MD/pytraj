#!bin/sh

# create this file to hide output
# python setup.py install --amber-release


if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    python setup.py build_ext -i --disable-openmp
    pip install .
else
    if [ "$TEST_SETUP" = 'true' ]; then
        echo "TEST_SETUP"
    else
        python setup.py build_ext -i
        pip install .
    fi
fi
