from ctypes.util import find_library
import ctypes
import os
import sys
from glob import glob

def find_libcpptraj():
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    libcpptraj_path_list = []

    for path in paths:
        path = path.strip()
        fnamelist = glob(os.path.join(path, "libcpptraj*"))
        for fname in fnamelist: 
            if os.path.isfile(fname):
                libcpptraj_path_list.append(fname)

    if not libcpptraj_path_list:
        raise ImportError('can not find libcpptraj. '
                          'make sure to update your LD_LIBRARY_PATH')
    else:
        return libcpptraj_path_list
