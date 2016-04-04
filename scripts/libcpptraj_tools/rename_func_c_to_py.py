import sys
from .util import func_c_to_py

fname = sys.argv[1]
try:
    fname2 = sys.argv[2]
except:
    fname2 = "tmp.pyx"

with open(fname, 'r') as fh, open(fname2, 'w') as fhwrite:
    lines = fh.readlines()
    func_c_to_py(lines)
    fhwrite.writelines(lines)
