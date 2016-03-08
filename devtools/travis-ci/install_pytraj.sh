#!bin/sh

# create this file to hide output
# python setup.py install --amber-release

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    python setup.py build --disable-openmp

    cd cpptraj
    export CPPTRAJHOME=`pwd`
    cd ..

    echo "update @rpath for libcpptraj"
    for x in $(find build/ -name '*.so'); do
        name="libcpptraj.dylib"
        install_name_tool -change "@rpath/$name" "$CPPTRAJHOME/lib/$name" $x
    done
    echo "done building"
fi

echo "install pytraj"
python setup.py install
