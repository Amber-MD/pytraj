#!/usr/bin/env python
'''Aim: 
    I need to write *.pxd header files for cython based on Cpptraj *.h header file.
    this script will print list of files that have not been written'''

from glob import glob

cpptrajsrc = "/mnt/raidc/haichit/AMBER14_official.naga84.forPythonTest/AmberTools/src/cpptraj/src/"

list = glob(cpptrajsrc + "*.h")
mypxdlist = glob("*.pxd")

for headercpp in glob(cpptrajsrc + "*.h"):
    header = headercpp.split("/")[-1]
    root = header.split(".")[0]
    tmppxd = root + ".pxd"
    if not tmppxd in mypxdlist:
        print root 
