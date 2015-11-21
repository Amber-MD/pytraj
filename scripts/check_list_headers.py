#!/usr/bin/env python
'''Aim:
    I need to write *.pxd header files for cython based on Cpptraj *.h header file.
    this script will print list of files that have not been written'''


import os
from glob import glob
from itertools import chain

pyx_include_dirs = [
    directory for directory, dirs, files in os.walk('pytraj')
    if '__init__.pyx' in files or '__init__.pxd' in files
    or '__init__.py' in files
]

pyx_include_patterns = [
    p + '/*.pyx' for p in pyx_include_dirs]

cpptrajsrc = os.environ['CPPTRAJHOME'] + "/src/"
print(cpptrajsrc)

_mypyxlist = [glob(x) for x in pyx_include_patterns]
mypyxlist = list(chain.from_iterable(_mypyxlist))
long_pyx_list = ".".join(mypyxlist)

for headercpp in glob(cpptrajsrc + "*.h"):
    header = headercpp.split("/")[-1]
    root = header.split(".")[0]
    tmppyx = root + ".pyx"
    #print (tmppyx)
    if not tmppyx in long_pyx_list:
        print(root)
        # pass
