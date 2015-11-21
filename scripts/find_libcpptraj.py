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

    # make a copy to avoid looping
    tmp_path = paths[:]

    for p in tmp_path:
        if 'anaconda' in p:
            paths.append(os.path.join(p.split("anaconda")[0], "anaconda", "lib"))
            paths.append(os.path.join(p.split("anaconda")[0], "anaconda2", "lib"))
            paths.append(os.path.join(p.split("anaconda")[0], "anaconda3", "lib"))

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
    print(find_lib('cpptraj', True))
