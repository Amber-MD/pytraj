from glob import glob
import sys
import os

keyword = sys.argv[1]

pyfiles = glob("*.pxd")

for pyx in pyfiles:
    with open(pyx, 'r') as fh:
        txt = fh.read()
        if keyword in txt:
            cppfile = pyx.split(".")[0] + ".cpp"
            os.remove(cppfile)
