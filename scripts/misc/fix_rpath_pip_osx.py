#!/usr/bin/env python

import os, sys
import shutil
import subprocess
from glob import glob
from conda_build import post
from conda_build.os_utils import macho
from auditwheel import wheeltools

# 1. Install wheel: pip wheel .
# 2. Use this script to fix wheel
# 3. Double-check: wheel unpack your_new.whl

# Note: your_new.whl will be in ./wheelhouse folder
version = '.'.join(str(i) for i in sys.version_info[:2])

LIBCPPTRAJ_RPATH = '@rpath/python{}/site-packages/pytraj/lib/libcpptraj.dylib'.format(version)

def copy_libcpptraj_to_pytraj_lib(libcpptraj):
    try:
        os.mkdir('pytraj/lib')
    except OSError:
        pass

    shutil.copy(libcpptraj, 'pytraj/lib')
    os.system('install_name_tool -id {} pytraj/lib/libcpptraj.dylib'.format(LIBCPPTRAJ_RPATH))


def main(pkg_name, whl_name, libcpptraj):
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
            os.system('ls pytraj')
            copy_libcpptraj_to_pytraj_lib(libcpptraj)
            for root, dirs, files in os.walk(pkg_name):
                for fn in (root + '/' +  _ for _ in files):
                    if fn.endswith('.so'):
                        libcpptraj_dir = ''
                        for lib in macho.get_dylibs(fn):
                            if 'libcpptraj' in lib:
                                libcpptraj_dir = lib
                                break
                        subprocess.check_call(['install_name_tool', '-change', 
                                               libcpptraj_dir,
                                               LIBCPPTRAJ_RPATH,
                                               fn])
                        os.system('otool -L {}'.format(fn))
if __name__ == '__main__':
    pkg_name = 'pytraj'
    whl_name = sys.argv[1]
    cpptrajhome = os.getenv('CPPTRAJHOME', os.path.abspath('../cpptraj'))
    libcpptraj = cpptrajhome + '/lib/libcpptraj.dylib'
    main(pkg_name, whl_name, libcpptraj)
