#!/usr/bin/env python

import os, sys
import shutil
import subprocess
from glob import glob
from auditwheel import wheeltools

# 1. Install wheel: pip wheel .
# 2. Use this script to fix wheel
# 3. Double-check: wheel unpack your_new.whl

# Note: your_new.whl will be in ./wheelhouse folder
def get_dylibs(fn):
    output = subprocess.check_output([
        'otool',
        '-L',
        fn
    ]).decode()
    return [line.split()[0] for line in output.split('\n') if line]

def copy_libcpptraj_to_pytraj_lib(libcpptraj, python_version):
    pytraj_lib = 'pytraj/lib'
    new_libcpptraj = pytraj_lib + '/libcpptraj.dylib'
    LIBCPPTRAJ_RPATH = '@rpath/python{}/site-packages/pytraj/lib/libcpptraj.dylib'.format(python_version)
    try:
        os.mkdir(pytraj_lib)
    except OSError:
        pass
    shutil.copy(libcpptraj, pytraj_lib)
    os.system('install_name_tool -id {} {}'.format(LIBCPPTRAJ_RPATH, new_libcpptraj))
    netcdf_libs = [lib for lib in get_dylibs(new_libcpptraj) if 'libnetcdf.7.dylib' in lib]

    if netcdf_libs:
        fn  = netcdf_libs[0]
        subprocess.check_call(['install_name_tool', '-change', 
                               fn,
                               '@rpath/libnetcdf.7.dylib',
                               new_libcpptraj])
    return LIBCPPTRAJ_RPATH


def main(pkg_name, whl_name, libcpptraj, python_version):
    try:
        os.mkdir('wheelhouse')
    except OSError:
        pass
    with wheeltools.InWheel(whl_name, out_wheel='wheelhouse/{}'.format(whl_name)):
        LIBCPPTRAJ_RPATH = copy_libcpptraj_to_pytraj_lib(libcpptraj, python_version)
        for root, dirs, files in os.walk(pkg_name):
            for fn in (root + '/' +  _ for _ in files):
                if fn.endswith('.so'):
                    libcpptraj_dir = ''
                    for lib in get_dylibs(fn):
                        if 'libcpptraj' in lib:
                            libcpptraj_dir = lib
                            break
                    subprocess.check_call(['install_name_tool', '-change', 
                                           libcpptraj_dir,
                                           LIBCPPTRAJ_RPATH,
                                           fn])
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Fix osx wheel')
    parser.add_argument('whl_name')
    parser.add_argument('--py', default='3.5', help='Python version')
    args = parser.parse_args()
    pkg_name = 'pytraj'
    whl_name = args.whl_name
    python_version = args.py
    cpptrajhome = os.getenv('CPPTRAJHOME', os.path.abspath('../cpptraj'))
    libcpptraj = cpptrajhome + '/lib/libcpptraj.dylib'
    main(pkg_name, whl_name, libcpptraj, python_version)
