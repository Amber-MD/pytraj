#!/usr/bin/env python

from __future__ import print_function
from subprocess import call
import sys
import os
from data.find_dupes import find_dupes

PY3 = sys.version_info[0] == 3
# for testing
i_duplicated = find_dupes(".") + find_dupes("./data")
print('number of duplicated filename (for mac) = ', i_duplicated)

with open("log", 'w') as log_file:
    with open("output.txt", 'w') as file_out:
        # get all the files starting with 'test_' and having "import unittest"
        call(['python', './get_unittest_files.py'])
        if PY3:
            with open("./TestListTravis.sh", 'r') as fh0:
                txt = fh0.read()
                txt = txt.replace("python", "python3")
            with open("./TestListTravis_py3.sh", 'w') as fh1:
                fh1.write(txt)

        # run tests
        if PY3:
            call(['sh', './/TestListTravis_py3.sh'],
                 stdout=file_out,
                 stderr=log_file)
        else:
            call(['sh', './/TestListTravis.sh'],
                 stdout=file_out,
                 stderr=log_file)

with open("log", 'r') as log_file, open("log2.sh", 'w') as log2:
    i_fails = 0  # only count files failed the assert
    i_seg = 0  # count Segmentation too
    for line in log_file.readlines():
        if 'File "./test_' in line:
            test = line.split()[1].replace(",", " ")  # "test_"
            log2.write("python " + test + "\n")
            i_fails += 1
        if 'Segmentation' in line:
            print(line)
            i_seg += 1

os.system("sh log2.sh")
print("%s FAILs (assert)" % i_fails)
print("no-passing-files: \n")
print("got %s segmentations faults" % i_seg)
os.system("cat log2.sh")

if i_fails >= 1 or i_seg >= 1 or i_duplicated >= 1:
    sys.exit(1)
