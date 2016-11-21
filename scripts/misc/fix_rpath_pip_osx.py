#!/usr/bin/env python

import os, sys
from conda_build import post
from glob import glob
from auditwheel import wheeltools

# 1. Install wheel: pip wheel .
# 2. Use this script to fix wheel
# 3. Double-check: wheel unpack your_new.whl

# Note: your_new.whl will be in ./wheelhouse folder

pkg_name = 'pytraj'
whl_name = sys.argv[1]

def copy_libcpptraj_to_pytraj_lib(libcpptraj):
    try:
        os.mkdir('pytraj/lib')
    except OSError:
        pass

    shutil.copy('{}/libcpptraj.dylib'.format(libcpptraj), 'pytraj/lib')
    os.system('install_name_tool -id @rpath/libcpptraj.dylib pytraj/lib/libcpptraj.dylib')


def main():
    if not whl_name.endswith('.whl'):
        pkg_name = whl_name
        for fn in (glob('{}/*so'.format(pkg_name)) + glob('{}/*/*so'.format(pkg_name))):
            post.mk_relative_osx(fn, pkg_name)
    else:
        try:
            os.mkdir('wheelhouse')
        except OSError:
            pass
        with wheeltools.InWheel(whl_name, out_wheel='wheelhouse/{}'.format(whl_name)):
            for fn in (glob('{}/*so'.format(pkg_name)) + glob('{}/*/*so'.format(pkg_name))):
                post.mk_relative_osx(fn, pkg_name)

if __name__ == '__main__':
    main()
