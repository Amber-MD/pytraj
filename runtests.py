#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import time

my_script = sys.argv[0]

try:
    verbose = sys.argv[1] in ['-verbose', 'verbose', '-v']
    need_help = sys.argv[1] in ['help', '-help', '--help']
    do_simple_test = sys.argv[1] in ['simple', 'minimal', '-simple',
                                     '-minimal', 'sim']
except:
    verbose = False
    need_help = False
    do_simple_test = False

if need_help:
    print("example")
    print("    verbose run: python %s -verbose" % my_script)
    print("    quiet run: python %s" % my_script)
    print("    short testing: python %s -simple" % my_script)
    sys.exit(0)

print("start testing. Go to ./tests folder")
os.chdir("./tests/")

if do_simple_test:
    os.system("python ./run_simple_test.py")
    sys.exit(0)

if verbose:
    os.system("python -m unittest")
else:
    os.system("python ./run_all_and_find_fails.py")

print("end testing")
