#!/usr/bin/env python
import sys
from six import string_types
"""replace old key word by new key word in a file and save it to ./tmp/"""
# for fh in `grep six_2 *py | sed "s/:/ /" | awk '{print $1}'`; do
# python ../scripts/replace_keyword.py six_2 compat $fh; done

oldword, newword = sys.argv[1:3]
oldfiles = sys.argv[3:]


if isinstance(oldfiles, string_types):
    oldfiles = [oldfiles, ]
else:
    oldfiles = oldfiles

for oldfile in oldfiles:
    with open(oldfile, 'r') as fh:
        txt = fh.read()

    if oldword in txt:
        txt = txt.replace(oldword, newword)

        with open("./tmp/" + oldfile, 'w') as fh2:
            fh2.write(txt)
    else:
        print("no word to replace")
