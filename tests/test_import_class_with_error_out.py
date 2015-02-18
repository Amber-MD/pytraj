#!/usr/bin/env python

import unittest
import sys
import os

modlist = []
with open("../pyxlist.txt") as pylist:
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
    #print "test import pytraj.%s" % mod
    classname = ".".join(['pytraj', mod])
    if "/" in classname:
        classname = classname.replace("/", ".")
    __import__(classname)
