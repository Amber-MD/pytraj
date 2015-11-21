#!/usr/bin/env python

from glob import glob
import os

lines = []
base_line = "import unittest"

# make a list of test cases
testlist = glob("test_*.py")
for folder in ['cluster',
               'test_frames',
               'test_trajs',
               'test_topology',
               'issues',
               'test_memory_usage',
               'test_iterators',
               'test_load', ]:
    testlist += glob(os.path.join(folder, 'test_*.py'))

for pyfile in testlist:
    with open(pyfile, 'r') as fh:
        txt = fh.read()
        # only test if file having keyword `import unittest`
        # but exclude ones with `#import unittest`
        # or one with "import unittest # pragma no_test"
        if base_line in txt:
            if not base_line + " # pragma no_test" in txt and not "#" + base_line in txt:
                line0 = "echo ./%s \n" % pyfile
                line = "python ./%s \n" % pyfile
                # print(pyfile)
                lines.append(line0)
                lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
