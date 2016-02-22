# this script was adapted from pythran package
# https://github.com/serge-sans-paille/pythran
# BSD-3 license
# https://github.com/serge-sans-paille/pythran/blob/master/LICENSE

from __future__ import print_function
import os
import sys
from ctypes.util import find_library
from subprocess import check_output

cxx = os.environ.get('CXX', 'c++')


def get_openmp_flag():
    """ Find OpenMP library and try to load if using ctype interface. """
    # find_library() does not search automatically LD_LIBRARY_PATH
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    for gomp in ('libgomp.so', 'libgomp.dylib'):
        cmd = [cxx, '-print-file-name=' + gomp]
        # the subprocess can fail in various ways
        # in that case just give up that path
        try:
            path = os.path.dirname(check_output(cmd).strip())
            if path:
                paths.append(path)
        except OSError:
            pass

    # Try to load find libgomp shared library using loader search dirs
    libgomp_path = find_library("gomp")

    # Try to use custom paths if lookup failed
    for path in paths:
        if libgomp_path:
            break
        path = path.strip()
        if os.path.isdir(path):
            libgomp_path = find_library(os.path.join(path, "libgomp"))

    if not libgomp_path:
        return ''
    else:
        return '-openmp'

if __name__ == '__main__':
    print(get_openmp_flag())
