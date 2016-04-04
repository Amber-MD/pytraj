import os
import sys
from codegen_enum import create_enum_of_dict

try:
    mode = sys.argv[1]
except:
    raise ValueError("must specify keyword for sys.argv[1]")

flist = ['MetaData',
         'DataSet',
         ]

for fname in flist:
    fname = fname + ".h"
    cpptrajsrc = os.environ['CPPTRAJHOME'] + "/src/"
    fname_full = cpptrajsrc + fname

    if os.path.exists(fname_full):
        create_enum_of_dict(fname, mode=mode)
