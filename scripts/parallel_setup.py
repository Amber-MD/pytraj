#!/usr/bin/env python

'''python parallel_setup.py

use 4 cores
Do not use this script for now, got
*** Error in `python': double free or corruption (out): 0x0000000002678b20 ***
'''

import os
from multiprocessing import Pool


def func(rank):
    # openmp
    os.system('python setup.py build_ext --rank=%s -i' % rank)

    # double free memory
    #os.system('python setup.py build --rank=%s' % rank)

# chose n_cpus=6 to be compromised with n_pyx_files=19
n_cpus = 6

Pool(n_cpus).map(func, range(n_cpus))

os.system('python setup.py install openmp')
