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
def get_so_files(folder):
    # all python extensions needed to be fixed
    level = 6
    all_so_files = []
    cmd = "find {} -type f -name '*.so'".format(folder)
    output = subprocess.check_output(cmd, shell=True)
    return [fn for fn in output.decode().split('\n') if fn]


def get_dylibs(fn):
    output = subprocess.check_output([
        'otool',
        '-L',
        fn
    ]).decode()
    return [line.split()[0] for line in output.split('\n') if line]


def copy_and_update_libs(pytraj_dir, libcpptraj, python_version):
    """ Copy libcpptraj.dylib and libnetcdf.7.dylib to pytraj/lib folder
    Update their ids to @rpath/{libcpptraj.dylib, libnetcdf.7.dylib}

    Also add loader_path to .so files in pytraj folder.
    """
    pytraj_lib = os.path.join(pytraj_dir, 'lib')
    new_libcpptraj = pytraj_lib + '/libcpptraj.dylib'

    try:
        os.mkdir(pytraj_lib)
    except OSError:
        pass
    shutil.copy(libcpptraj, pytraj_lib)
    netcdf_lib = [lib for lib in get_dylibs(new_libcpptraj) if 'libnetcdf.7.dylib' in lib][0]
    shutil.copy(netcdf_lib, pytraj_lib)
    os.system('install_name_tool -id {} {}'.format('@rpath/libcpptraj.dylib', new_libcpptraj))
    os.system('install_name_tool -id {} {}'.format('@rpath/libnetcdf.7.dylib', netcdf_lib))


def main(pkg_name, whl_name, libcpptraj, python_version):
    try:
        os.mkdir('wheelhouse')
    except OSError:
        pass
    with wheeltools.InWheel(whl_name, out_wheel='wheelhouse/{}'.format(whl_name)):
        LIBCPPTRAJ_RPATH = copy_and_update_libs(libcpptraj, python_version)
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
