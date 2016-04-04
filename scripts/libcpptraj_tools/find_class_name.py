# (c) 2014 Hai Nguyen
#!/usr/bin/env python
# python ./find_class_name.py filename

import os
import sys
import re
from glob import glob

CPPTRAJSRC = os.environ['CPPTRAJHOME'] + "/src/"


def find_class(src=CPPTRAJSRC):
    #src = os.environ['AMBERHOME'] + '/AmberTools/src/cpptraj/src/'
    p = re.compile(r'#include "(.+?).h"')
    classlist = []

    for fname in glob(src + "*.h"):
        fnshort = fname.split("/")[-1]
        fh = open(fname, 'r')
        for line in fh.readlines():
            if line.startswith("class"):
                classname = line.split()[1].split(":")[0].split(";")[0]
                classlist.append(classname)
        fh.close()
    return list(set(classlist))

if __name__ == '__main__':
    classlist = find_class()
    for x in classlist:
        if (not x.startswith("Action")) and (not x.startswith("Analysis")):
            print(x)
