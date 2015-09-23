#!/usr/bin/env python

import os
import sys


def find_dupes(my_dir='.'):
    myfiles = set(os.listdir(my_dir))
    myfiles_lower = set()

    i_duplicated = 0
    for myfile in myfiles:
        lower = myfile.lower()
        if lower in myfiles_lower:
            i_duplicated += 1
            print('\033[91m%s\033[0m differs only by case! Change its name.' %
                  myfile)
        else:
            myfiles_lower.add(lower)
    return i_duplicated


if __name__ == "__main__":
    n = find_dupes()
    print(n)
