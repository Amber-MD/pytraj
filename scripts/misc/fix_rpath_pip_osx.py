#!/usr/bin/env python

import os, sys
import shutil
import subprocess
from glob import glob
from auditwheel import wheeltools

# 1. Install wheel: pip wheel .
# 2. Use this script to fix wheel
# 3. Double-check: wheel unpack your_new.whl


def get_dylibs(fn):
    output = subprocess.check_output([
        'otool',
        '-L',
        fn
    ]).decode()
    return [line.split()[0] for line in output.split('\n') if line]


def copy_and_update_libs(pytraj_dir, libcpptraj):
    """ Copy libcpptraj.dylib and libnetcdf.7.dylib to pytraj/lib folder
    Update their ids to @rpath/{libcpptraj.dylib, libnetcdf.7.dylib}

    Also add loader_path to .so files in pytraj folder.
    """
    pytraj_lib = os.path.join(pytraj_dir, 'pytraj_3rd_party')
    new_libcpptraj = pytraj_lib + '/libcpptraj.dylib'
    new_libnetcdf = pytraj_lib + '/libnetcdf.7.dylib'

    try:
        os.mkdir(pytraj_lib)
    except OSError:
        pass
    shutil.copy(libcpptraj, new_libcpptraj)
    netcdf_lib = [lib for lib in get_dylibs(libcpptraj) if 'libnetcdf.7.dylib' in lib][0]
    print('netcdf_lib', netcdf_lib)
    shutil.copy(netcdf_lib, new_libnetcdf)
    os.system('install_name_tool -id {} {}'.format('@rpath/pytraj_3rd_party/libcpptraj.dylib', new_libcpptraj))
    os.system('sudo install_name_tool -id {} {}'.format('@rpath/pytraj_3rd_party/libnetcdf.7.dylib', new_libnetcdf))


def pack(pkg_name, whl_name, libcpptraj):
    print('libcpptraj', libcpptraj)
    try:
        os.mkdir('wheelhouse')
    except OSError:
        pass
    with wheeltools.InWheel(whl_name, out_wheel='wheelhouse/{}'.format(whl_name)):
        pkg_dir = os.path.abspath(pkg_name)
        copy_and_update_libs(pkg_dir, libcpptraj)
        rpath_libcpptraj = '@rpath/pytraj_3rd_party/libcpptraj.dylib'

        lib_dir = os.path.join(pkg_name, 'lib')

        for root, dirs, files in os.walk(pkg_name):
            for fn in (os.path.join(root,  _) for _ in files):
                if fn.endswith('.so'):
                    libcpptraj_dir = ''
                    for lib in get_dylibs(fn):
                        if 'libcpptraj' in lib:
                            libcpptraj_dir = lib
                            break
                    relpath = os.path.relpath(pkg_dir, os.path.abspath(os.path.dirname(fn)))
                    print('relpath', relpath)
                    subprocess.check_call(['install_name_tool', '-change', 
                                           libcpptraj_dir,
                                           rpath_libcpptraj,
                                           fn])
                    loader_path = '@loader_path/{}/'.format(relpath)
                    print('loader_path', loader_path)
                    subprocess.check_call(['install_name_tool', '-add_rpath', 
                                           loader_path,
                                           fn])

def main():
    import argparse
    parser = argparse.ArgumentParser('Fix osx wheel')
    parser.add_argument('whl_name')
    parser.add_argument('--libcpptraj')
    args = parser.parse_args()
    pkg_name = 'pytraj'
    pack(pkg_name, args.whl_name, os.path.abspath(args.libcpptraj))


if __name__ == '__main__':
    main()
