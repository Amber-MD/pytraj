from __future__ import print_function

import os
from glob import glob

def find_lib(libname, unique=False):
    """return a list of all library files"""
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    lib_path_list = []
    key = "lib" + libname + "*"

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
    print (find_lib('cpptraj', True))
