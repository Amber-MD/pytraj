#!/bin/sh 
pyxfile=$1

cat >setup_single_file.py <<EOF
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

pyxfile = "$pyxfile"
ext = pyxfile.split(".")[0]

setup(
      ext_modules = cythonize([
          Extension(ext, ["$pyxfile",],
          ])))
          
EOF

python ./setup_single_file.py build_ext -i
