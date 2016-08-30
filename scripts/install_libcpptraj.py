#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import subprocess
sys.path.append('scripts')
from check_openmp import get_openmp_flag
from find_lib import find_lib


# DEFAULT_MAC_BUILD = '-shared -macAccelerate --with-fftw3=/usr/local --with-netcdf=/usr/local -noarpack'
DEFAULT_MAC_BUILD = '-shared -macAccelerate -noarpack'

DEFAULT_MAC_CCOMPILER = 'clang'
DEFAULT_MAC_CXXCOMPILER = 'clang++'

CPPTRAJ_CXX = ' '
IS_OSX = (sys.platform == 'darwin')

if sys.platform.startswith('darwin'):
    # Hack to fix conda
    import distutils.sysconfig as sc
    compiler = sc.get_config_var('CC') or ''
    if compiler == 'gcc':
        # This must be a conda install that uses a a pure compiler name (rather
        # than an absolute path).
        DEFAULT_MAC_CCOMPILER = "/usr/bin/gcc"
        DEFAULT_MAC_CXXCOMPILER = "/usr/bin/g++"

        # Warning: dirty hack to use libstdc++
        CPPTRAJ_CXX = '-stdlib=libstdc++'

def add_CPPTRAJ_CXX_to_config(fn, CPPTRAJ_CXX=CPPTRAJ_CXX):
    with open(fn, 'r') as fconfig:
        lines = fconfig.readlines()

    with open('tmp.h', 'w') as fh:
        for idx, line in enumerate(lines):
            if line.startswith('CXX='):
                print('line', line)
                lines[idx] = ' '.join((line.strip(), CPPTRAJ_CXX, '\n'))
        fh.write(''.join(lines))
    subprocess.check_call('mv tmp.h {}'.format(fn), shell=True)


def get_compiler_and_build_flag():
    args = parse_args()

    openmp_flag = '-openmp' if args.openmp else ''
    if openmp_flag:
        assert get_openmp_flag(), 'your system must support openmp'

    install_type = args.install_type
    debug = '-debug' if args.debug else ''

    try:
        import numpy as np
        has_numpy = True
    except ImportError:
        has_numpy = False

    # better name for CPPTRAJ_COMPILER_OPTION?
    # cpptraj: ./configure gnu
    # e.g: CPPTRAJ_COMPILER_OPTION=gnu python ./scripts/install_libcpptraj.py
    cpptraj_compiler_option = os.environ.get('CPPTRAJ_COMPILER_OPTION', 'gnu')  # intel | pgi | clang | cray?
    amberhome = os.environ.get('AMBERHOME', '')
    amberlib = '-amberlib' if amberhome and args.amberlib else ''

    if has_numpy and find_lib('openblas'):
        prefix = sys.prefix
        # likely having openblas?
        build_flag_ = ('--with-netcdf={prefix} --with-blas={prefix} '
                       '--with-bzlib={prefix} --with-zlib={prefix} '
                       '-openblas -noarpack'.format(prefix=prefix))
    elif has_numpy:
        try:
            blas_prefix = np.__config__.blas_opt_info['library_dirs'][0].strip('lib')
            lapack_prefix = np.__config__.lapack_opt_info['library_dirs'][0].strip('lib')
            build_flag_ = ('-noarpack --with-blas={blas_prefix} --with-lapack={lapack_prefix}'
                .format(blas_prefix=blas_prefix, lapack_prefix=lapack_prefix))
        except (KeyError, IndexError):
            build_flag_ = '-noarpack'
    else:
        # user gets lucky?
        build_flag_ = '-noarpack'

    if IS_OSX:
        build_flag = ' '.join((DEFAULT_MAC_BUILD, amberlib))
    else:
        build_flag = ' '.join(('-shared', build_flag_, amberlib, openmp_flag))

    build_flag = ' '.join((build_flag, debug, '-nosanderlib'))

    if install_type == 'github':
        print('install libcpptraj from github')
        subprocess.check_call('git clone https://github.com/Amber-MD/cpptraj'.split())
    else:
        print('install libcpptraj from current ./cpptraj folder')
    return cpptraj_compiler_option, build_flag

def fix_rpath():
    if IS_OSX:
        subprocess.check_call('(cd lib && install_name_tool -id `pwd`/libcpptraj.dylib libcpptraj.dylib)'.split())
    else:
        pass

def install_libcpptraj(cpptraj_compiler_option, build_flag):
    cwd = os.getcwd()
    try:
        os.chdir('./cpptraj')
    except OSError:
        raise OSError('please try python scripts/install_cpptraj.py github')
    os.environ['CPPTRAJHOME'] = os.getcwd()

    try:
        os.mkdir('./lib')
    except OSError:
        pass


    if cpptraj_compiler_option == 'clang' and sys.platform == 'darwin':
        cxx_overwrite = 'CXX="clang++ -stdlib=libstdc++"'
    else:
        cxx_overwrite = ''

    cm = 'bash configure {build_flag} {compiler} {cxx_overwrite}'.format(
            build_flag=build_flag, compiler=cpptraj_compiler_option, cxx_overwrite=cxx_overwrite)

    print('build command: ', cm)
    subprocess.check_call(cm.split())

    if IS_OSX:
        add_CPPTRAJ_CXX_to_config('config.h', CPPTRAJ_CXX)

    subprocess.check_call('make libcpptraj -j8'.split())
    fix_rpath()
    os.chdir(cwd)

    print("make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
          "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
          "then 'python ./setup.py install'")

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='install options')
    parser.add_argument('-openmp', action='store_true', help='use openmp')
    parser.add_argument('-amberlib', action='store_true', help='use use amberlib if $AMBERHOME is set')
    parser.add_argument('-debug', action='store_true', help='debug')
    parser.add_argument('install_type', default='', nargs='?', help='install_type in amber')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    cpptraj_compiler_option, build_flag = get_compiler_and_build_flag()
    install_libcpptraj(cpptraj_compiler_option, build_flag)
