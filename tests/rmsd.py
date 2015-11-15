#!/usr/bin/env python

import sys
from numpy import sum, sqrt, loadtxt


def rmsd(a1, a2):
    '''rmsd of a1 and a2
    '''
    alen = len(a1)
    if alen != len(a2):
        print("two array must have the same length")
        sys.exit()
    tmp = sum((a1 - a2)**2)
    return sqrt(tmp / alen)


if __name__ == '__main__':
    a1, a2 = sys.argv[1:3]
    a1_arr = loadtxt(a1)
    a2_arr = loadtxt(a2)
    print(rmsd(a1_arr, a2_arr))
