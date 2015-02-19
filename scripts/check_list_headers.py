#!/usr/bin/env python
'''Aim: 
    I need to write *.pxd header files for cython based on Cpptraj *.h header file.
    this script will print list of files that have not been written'''


import os
from glob import glob

cpptrajsrc = os.environ['CPPTRAJHOME'] + "/src"
print (cpptrajsrc)

list = glob(cpptrajsrc + "*.h")
mypxdlist = glob("*.pxd")

for headercpp in glob(cpptrajsrc + "*.h"):
    header = headercpp.split("/")[-1]
    root = header.split(".")[0]
    tmppxd = root + ".pxd"
    if not tmppxd in mypxdlist:
        print (root) 
