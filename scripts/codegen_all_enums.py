import os
import sys
from codegen_enum import create_enum_of_dict
from find_class_name import find_class

try:
    mode = sys.argv[1]
except:
    raise ValueError("must specify keyword for sys.argv[1]")

for fname in find_class():
    fname = fname + ".h"
    #cpptrajsrc = os.environ['AMBERHOME'] + "AmberTools/src/cpptraj/src/"
    cpptrajsrc = os.environ['CPPTRAJHOME'] + "/src/"
    fname_full = cpptrajsrc + fname
    if os.path.exists(fname_full):
        create_enum_of_dict(fname, mode=mode)
