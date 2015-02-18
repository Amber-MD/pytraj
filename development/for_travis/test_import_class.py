#!/usr/bin/env python

import sys
import os

modlist = []
PYCPPTRAJ_HOME = os.getcwd()
with open(PYCPPTRAJ_HOME + "/PYXLIST.txt") as pylist:
    lines = pylist.readlines()
    for line in lines:
        if "#" not in line:
            modlist.append(line.split()[0])

try:
    # if exist sys.argv[1]: test only this mod name
    modlist = [sys.argv[1],]
except:
    # else: use modlist just created above
    pass

for mod in modlist:
    #print "test import pycpptraj.%s" % mod
    classname = ".".join(['pycpptraj', mod])
    try:
        __import__(classname)
    except ImportError:
        print "Error mod = %s" % classname

print "end test importing"
