from glob import glob
from util import _to_lower_case
import re

need2fixed = []
word = "object has no attribute"
with open("log", 'r') as fh_log:
    for line in fh_log.readlines():
        if word in line:
            tmp = line.split(word)[-1]
            tmp = tmp.replace("'", "")
            need2fixed.append(tmp.split()[0])
print(need2fixed)
newlist = [_to_lower_case(x) for x in need2fixed]
print(newlist)

for testfile in glob("test*.py"):
    found_word = False
    with open(testfile, 'r') as fhtest:
        lines = fhtest.readlines()
        for i, line in enumerate(lines):
            for oldword, newword in zip(need2fixed, newlist):
                tmp = ".%s" % oldword
                pytmp = "pytraj" + tmp
                if tmp in line and pytmp not in line:
                    lines[i] = line.replace(tmp, "." + newword)
                    found_word = True
    if found_word:
        with open(testfile + "_", 'w') as fh:
            fh.writelines(lines)
