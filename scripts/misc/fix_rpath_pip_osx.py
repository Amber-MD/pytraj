#!/usr/bin/env python

import os, sys
import shutil
import subprocess
from glob import glob
from auditwheel import wheeltools

# 1. Install wheel: pip wheel .
# 2. Use this script to fix wheel
# 3. Double-check: wheel unpack your_new.whl


def get_dylibs(fn, local_only=False):
    # local_only: /usr/local
    output = subprocess.check_output([
        'otool',
        '-L',
        fn
    ]).decode()
    lines = [line.split()[0] for line in output.split('\n') if line]
    if local_only:
        return [line for line in lines if '/usr/local' in line]
    else:
        return lines


def update_id(fn, sudo=False):
    # e.g: @rpath/libcpptraj.dylib
    basename = os.path.basename(fn)
    rpath = '@rpath/{}'.format(basename)
    cmd = 'install_name_tool -id {} {}'.format(rpath, fn)
    if sudo:
        cmd = 'sudo '  + cmd
    os.system(cmd)
    return rpath

def copy_and_update_libs(pytraj_dir, libcpptraj):
    """ Copy libcpptraj.dylib and libnetcdf.{version}.dylib to pytraj/lib folder
    Update their ids to @rpath/{libcpptraj.dylib, libnetcdf.{version}.dylib}

    Also add loader_path to .so files in pytraj folder.
    """
    pytraj_lib = os.path.join(pytraj_dir, 'pytraj_3rd_party')
    netcdf_lib = [lib for lib in get_dylibs(libcpptraj) if 'libnetcdf' in lib][0]
    netcdf_basename = os.path.basename(netcdf_lib)
    new_libcpptraj = pytraj_lib + '/libcpptraj.dylib'
    new_libnetcdf = os.path.join(pytraj_lib, netcdf_basename)

    try:
        os.mkdir(pytraj_lib)
    except OSError:
        pass
    shutil.copy(libcpptraj, new_libcpptraj)
    rpath_libnetcdf = '@rpath/{}'.format(netcdf_basename)
    print('netcdf_lib', netcdf_lib)
    shutil.copy(netcdf_lib, new_libnetcdf)
    os.system('install_name_tool -id {} {}'.format('@rpath/pytraj_3rd_party/libcpptraj.dylib', new_libcpptraj))
    os.system('install_name_tool -change {} {} {}'.format(netcdf_lib, rpath_libnetcdf,  new_libcpptraj))
    os.system('install_name_tool -add_rpath @loader_path/ {}'.format(new_libcpptraj))
    os.system('sudo install_name_tool -id {} {}'.format(rpath_libnetcdf, new_libnetcdf))

    for fn in get_dylibs(netcdf_lib, local_only=True):
        if netcdf_basename not in fn:
            new_fn = os.path.join(pytraj_lib, os.path.basename(fn))
            shutil.copy(fn, new_fn)
            rpath = update_id(new_fn, sudo=True)
            os.system('sudo install_name_tool -change {} {} {}'.format(fn, rpath, new_libnetcdf))

            for fn2 in get_dylibs(fn, local_only=True):
                rpath = '@rpath/{}'.format(os.path.basename(fn2))
                os.system('sudo install_name_tool -change {} {} {}'.format(fn2, rpath, new_fn))


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

        so_files = []
        for root, dirs, files in os.walk(pkg_name):
            for fn in (os.path.join(root,  _) for _ in files):
                if fn.endswith('.so'):
                    so_files.append(os.path.abspath(fn))
                    libcpptraj_dir = ''
                    for lib in get_dylibs(fn):
                        if 'libcpptraj' in lib:
                            libcpptraj_dir = lib
                            break
                    relpath = os.path.relpath(pkg_dir, os.path.abspath(os.path.dirname(fn)))
                    subprocess.check_call(['install_name_tool', '-change', 
                                           libcpptraj_dir,
                                           rpath_libcpptraj,
                                           fn])
                    loader_path = '@loader_path/{}/'.format(relpath)
                    subprocess.check_call(['install_name_tool', '-add_rpath', 
                                           loader_path,
                                           fn])
        print("rpath")
        for fn in glob(os.path.join(pkg_dir, 'pytraj_3rd_party/*dylib')):
            print(os.path.basename(fn))
            print("    ", get_dylibs(fn)[1:])
        print(get_dylibs(so_files[0])[1:])

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
