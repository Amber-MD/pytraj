from __future__ import print_function

import os
import sys
from glob import glob


def find_lib(libname, unique=False):
    """return a list of all library files"""
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    paths += os.environ.get('AMBERHOME', '').split(':')
    paths += os.environ.get('PYTHONPATH', '').split(':')
    paths += os.environ.get('CPPTRAJHOME', '').split(':')
    paths += os.environ.get('PATH', '').split(':')
    paths += os.environ.get('ANACONDAHOME', '').split(':')

    anconda_dir = os.environ.get('ANCONDAHOME', '')

    if anconda_dir:
        pattern = os.path.join(anconda_dir, 'envs', '*')
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
    #print(find_lib('cpptraj', True))
    #print(find_lib('libcpptraj', True))
    print(find_lib('libcpptraj.so', True))
