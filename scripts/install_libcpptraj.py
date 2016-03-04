#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import subprocess
from subprocess import CalledProcessError
sys.path.append('scripts')
from check_openmp import get_openmp_flag
from find_lib import find_lib


DEFAULT_MAC_BUILD = '-shared -macAccelerate --with-fftw3=/usr/local --with-netcdf=/usr/local -noarpack'
DEFAULT_MAC_COMPILER = 'clang'


def get_compiler_and_build_flag():
    try:
        sys.argv.remove('-openmp')
        openmp_flag = '-openmp'
        assert get_openmp_flag(), 'your system must support openmp'
    except ValueError:
        openmp_flag = ''
    
    try:
        install_type = sys.argv[1]
    except IndexError:
        install_type = ''
    
    try:
        import numpy as np
        has_numpy = True
    except ImportError:
        has_numpy = False

    compiler = os.environ.get('COMPILER', 'gnu')
    if sys.platform == 'darwin':
        compiler = DEFAULT_MAC_COMPILER
    amberhome = os.environ.get('AMBERHOME', '')
    amberlib = '-amberlib' if amberhome else ''
    
    if has_numpy and find_lib('openblas'):
        prefix = sys.prefix
        # likely having openblas?
        build_flag_ = '--with-netcdf={prefix} --with-blas={prefix} \
                       --with-bzlib={prefix} --with-zlib={prefix} \
                       -openblas -noarpack'.format(prefix=prefix)
    elif has_numpy:
        try:
            blas_prefix = np.__config__.blas_opt_info['library_dirs'][0].strip('lib')
            lapack_prefix = np.__config__.lapack_opt_info['library_dirs'][0].strip('lib')
            build_flag_ = '-noarpack --with-blas={blas_prefix} --with-lapack={lapack_prefix}'.format(blas_prefix=blas_prefix, lapack_prefix=lapack_prefix)
        except (KeyError, IndexError):
            build_flag_ = '-noarpack'
    else:
        # user gets lucky?
        build_flag_ = '-noarpack'
    
    if sys.platform == 'darwin':
        build_flag = DEFAULT_MAC_BUILD
    else:
        build_flag = ' '.join(('-shared', build_flag_, amberlib, openmp_flag))
    
    if install_type == 'github':
        print('install libcpptraj from github')
        os.system('git clone https://github.com/Amber-MD/cpptraj')
    else:
        print('install libcpptraj from current ./cpptraj folder')
    return compiler, build_flag

def install_libcpptraj(compiler, build_flag):
    cwd = os.getcwd()
    try:
        os.chdir('./cpptraj')
    except FileNotFoundError:
        raise FileNotFoundError('please try python scripts/install_cpptraj.py github')
    os.environ['CPPTRAJHOME'] = os.getcwd()
    
    try:
        os.mkdir('./lib')
    except FileExistsError:
        pass
    
    print('build_flag = ', build_flag)
    os.system('bash configure {build_flag} {compiler} || exit 1'.format(build_flag=build_flag, compiler=compiler))
    # os.system('make libcpptraj -j8 || exit 1')
    os.system('make libcpptraj || exit 1')
    os.chdir(cwd)
    print("make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
          "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
          "then 'python ./setup.py install'")


if __name__ == '__main__':
    compiler, build_flag = get_compiler_and_build_flag()
    install_libcpptraj(compiler, build_flag)
