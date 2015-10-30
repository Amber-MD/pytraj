# usr/bin/env python
from subprocess import call
import sys
import os
from glob import glob
from pytraj.compat import set

keyword = sys.argv[1]

keyword_flist = set(glob('*' + keyword + "*"))
all_test = set(glob("test_*.py"))

flist = keyword_flist & all_test

for fname in flist:
    print(fname)
    os.system("python %s" % fname)
