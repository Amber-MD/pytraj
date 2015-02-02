#!/bin/sh 
pyxfile=$1

cat >setup_single_file.py <<EOF
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

PXD_include = os.environ['PYCPPTRAJ_HOME'] + "pycpptraj"
PXD_include_action = os.environ['PYCPPTRAJ_HOME'] + "pycpptraj/Action/"
PXD_include_ana = os.environ['PYCPPTRAJ_HOME'] + "pycpptraj/Analysis/"
cpptraj_include = os.environ['AMBERHOME'] + "/AmberTools/src/cpptraj/src/"
amber_include =  os.environ['AMBERHOME'] + "/include/"
lib_dir = os.environ['AMBERHOME'] + "/lib/"
pyxfile = "$pyxfile"
ext = pyxfile.split(".")[0]

setup(
      ext_modules = cythonize([
          Extension(ext, ["$pyxfile",],
                    libraries=['cpptraj'],
                    library_dirs=[lib_dir,],
                    include_dirs=[PXD_include, PXD_include_ana, PXD_include_action, cpptraj_include, amber_include,])
          ])
     ) 
EOF

python ./setup_single_file.py build_ext -i
