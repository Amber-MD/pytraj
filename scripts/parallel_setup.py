#!/usr/bin/env python

'''python parallel_setup.py

use 4 cores
Do not use this script for now, got
*** Error in `python': double free or corruption (out): 0x0000000002678b20 ***
'''

import os
from multiprocessing import Pool

def func(rank):
    os.system('python setup.py build_ext --rank=%s' % rank)

Pool(4).map(func, range(4))

os.system('python setup.py install')
