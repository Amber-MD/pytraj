#!/bin/sh

isosx=`python -c "import sys; print(sys.platform.startswith('darwin'))"`
pyver=`python -c "import sys; print('.'.join(str(x) for x in sys.version_info[:2]))"`
pyver2=`python -c "import sys; print(''.join(str(x) for x in sys.version_info[:2]))"`

# FIXME: why hard code the version here?
version='2.0.4'

if [ "$isosx" = "True" ]; then
   whlfile=pytraj-$version-cp$pyver2-cp${pyver2}m-macosx_*_x86_64.whl
else
   if [ "$pyver2" = "27" ]; then
       whlfile=pytraj-$version-cp$pyver2-cp${pyver2}mu-manylinux1_x86_64.whl
   else
       whlfile=pytraj-$version-cp$pyver2-cp${pyver2}m-manylinux1_x86_64.whl
   fi
fi

wheel unpack `ls $RECIPE_DIR/../../../dist/wheelhouse/$whlfile`
cp -rf  pytraj-$version/pytraj* $PREFIX/lib/python$pyver/site-packages/
