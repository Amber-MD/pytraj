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
        # install_name_tool -change $name $CPPTRAJHOME/lib/$name $x
        # install_name_tool -change "@rpath/$name" "$CPPTRAJHOME/lib/$name" $x
        # install_name_tool -id "@rpath/$name" "$CPPTRAJHOME/lib/$name"
    done
    echo "done building"
    echo "install pytraj"
    python setup.py install --disable-openmp
    echo "copying libcpptraj"
    # sudo cp $CPPTRAJHOME/lib/libcpptraj.dylib /usr/local/lib/
else
    python setup.py install
fi

