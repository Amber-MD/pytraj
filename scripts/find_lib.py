'''
Example: python find_lib.py sander
'''

#!/usr/bin/env python

from __future__ import print_function

import os
import sys
from glob import glob
from itertools import chain


def find_lib(libname, unique=True):
    """return a list of all library files"""

    envlist = ['LD_LIBRARY_PATH', 'AMBERHOME', 'PYTHONPATH',
               'CPPTRAJHOME', 'PATH', 'ANCONDAHOME']
    paths = list(chain.from_iterable([os.environ.get(env_name, '').split(':') for env_name in envlist]))

    prefix_real_path = os.path.realpath(sys.prefix)
    if 'miniconda' in prefix_real_path:
        miniconda_dir = prefix_real_path
        paths.append(miniconda_dir + '/lib/')
    else:
        miniconda_dir = ''
        paths.append(prefix_real_path+ '/lib/')

    anconda_dir = os.environ.get('ANCONDAHOME', '')

    if anconda_dir:
        pattern = os.path.join(anconda_dir, 'envs', '*')
        paths += [os.path.join(dir, 'lib')  for dir in glob(pattern)]
    if miniconda_dir:
        pattern = os.path.join(miniconda_dir, 'envs', '*')
        paths += [os.path.join(dir, 'lib')  for dir in glob(pattern)]

    lib_path_list = []
    key = libname + '*' if libname.startswith('lib') else "lib" + libname + "*"

    for path in paths:
        path = path.strip()
        fnamelist = glob(os.path.join(path, key))
        for fname in fnamelist:
            if os.path.isfile(fname):
                lib_path_list.append(fname)

    if not lib_path_list:
        return None
    else:
        if unique:
            return set(lib_path_list)
        else:
            return lib_path_list


if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:
        sys.exit(1)

    print(find_lib(sys.argv[1]))
