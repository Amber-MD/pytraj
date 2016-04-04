import os
from util import Line_codegen
import sys

#cpp_src = os.environ['AMBERHOME'] + "AmberTools/src/cpptraj/src/"
cpp_src = "../cpptraj/src/"
fname = cpp_src + sys.argv[1]

with open(fname, 'r') as fh:
    lines = fh.readlines()
    for line in lines:
        if "typedef" in line:
            linegen = Line_codegen(line)
            linegen.remove_std_namespace()
            linegen.replace(";", "")
            linegen.replace("typedef", "ctypedef")
            linegen.replace_waka()
            print(linegen.myline)
