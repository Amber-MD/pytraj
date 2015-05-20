#!/usr/bin/env python

import os

myfiles = set(os.listdir('.'))
myfiles_lower = set()

for myfile in myfiles:
    lower = myfile.lower()
    if lower in myfiles_lower:
        print('\033[91m%s\033[0m differs only by case! Change its name.' % myfile)
    else:
        myfiles_lower.add(lower)
