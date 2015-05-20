#!/usr/bin/env python

import os
import sys

myfiles = set(os.listdir('.'))
myfiles_lower = set()

i_duplicated = 0
for myfile in myfiles:
    lower = myfile.lower()
    if lower in myfiles_lower:
        i_duplicated += 1
        print('\033[91m%s\033[0m differs only by case! Change its name.' % myfile)
    else:
        myfiles_lower.add(lower)

if i_duplicated > 0:
    sys.exit(1)
else:
    print ("no duplicated filenames")
