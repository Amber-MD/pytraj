#!/usr/bin/env  python

from subprocess import call
import sys
import os

PY3 = sys.version_info[0] == 3

with open("log", 'w') as log_file:
    with open("output.txt", 'w') as file_out:
        # get all the files starting with 'test_' and having "import unittest"
        call(['python', './get_py_files.py'])
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

with open("log", 'r') as log_file:
    i_fails = 0
    for line in log_file.readlines():
        if "File" in line:
            i_fails += 1
            print(line)

print("%s FAILs" % i_fails)
print("fail files: \n")
os.system("grep File log")

if i_fails >= 1:
    sys.exit(1)
