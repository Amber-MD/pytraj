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
IS_OSX = sys.platform.startswith('darwin')

if IS_OSX:
    # Hack to fix conda
    import distutils.sysconfig as sc
    compiler = sc.get_config_var('CC') or ''
    if compiler == 'gcc':
        # This must be a conda install that uses a a pure compiler name (rather
        # than an absolute path).
        DEFAULT_MAC_CCOMPILER = "/usr/bin/gcc"
        DEFAULT_MAC_CXXCOMPILER = "/usr/bin/g++"


def add_cpptraj_cxx_to_config(fn, cpptraj_cxx=CPPTRAJ_CXX):
    with open(fn, 'r') as fconfig:
        lines = fconfig.readlines()

    with open('tmp.h', 'w') as fh:
        for idx, line in enumerate(lines):
            if line.startswith('CXX='):
                print('line', line)
                lines[idx] = ' '.join((line.strip(), cpptraj_cxx, '\n'))
        fh.write(''.join(lines))
    subprocess.check_call('mv tmp.h {}'.format(fn), shell=True)

def ensure_gnu():
    # both cpptraj and pytraj give priority for CXX and CC environments
    # check them first.
    cc = os.getenv('CC', '')
    gcc_exe = cc if cc else 'gcc'
    out = subprocess.check_output([gcc_exe, '--version']).decode()
    if 'clang' in out:
        print('{} --version'.format(gcc_exe))
        print(out)
        print('{} here is actually clang compiler. Please export correct PATH for the real g++\n'.format(gcc_exe))
        print('Or export CXX and CC environments')
        print('e.g: ')
        print('    export CC=/usr/local/bin/gcc-5')
        print('    export CXX=/usr/local/bin/g++-5')
        sys.exit(1)

def get_compiler_and_build_flag():
    args = parse_args()

    openmp_flag = '-openmp' if args.openmp else ''
    if openmp_flag:
        # we turn off OSX openmp anyway, so do not need to check
        if not IS_OSX:
            assert get_openmp_flag(), 'your system must support openmp'

    install_type = args.install_type
    debug = '-debug' if args.debug else ''

    try:
        import numpy as np
        has_numpy = True
    except ImportError:
        has_numpy = False

    # better name for COMPILER?
    # cpptraj: ./configure gnu
    # e.g: COMPILER=gnu python ./scripts/install_libcpptraj.py
    os_dependent_default_compiler = 'clang' if IS_OSX else 'gnu'
    compiler_env = os.environ.get('COMPILER', '')  # intel | pgi | clang | cray?
    if compiler_env:
        print('libcpptraj: Using COMPILER env = ', compiler_env)
        compiler = compiler_env
    else:
        compiler = os_dependent_default_compiler
        cxx = os.getenv('CXX')
        if cxx:
            print('Using preset CXX', cxx)
            compiler = ''
    if compiler == 'gnu':
        ensure_gnu()
    amberhome = os.environ.get('AMBERHOME', '')
    amberlib = '-amberlib' if amberhome and args.amberlib else ''

    prefix = sys.prefix
    if has_numpy and find_lib('openblas'):
        # likely having openblas?
        build_flag_ = ('--with-netcdf={prefix} --with-blas={prefix} '
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

    zip_stuff = ' --with-bzlib={prefix} --with-zlib={prefix} '.format(prefix=prefix)
    build_flag_ += zip_stuff

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
    return compiler, build_flag

def install_libcpptraj(compiler, build_flag):
    '''

    Parameters
    ----------
    compiler : str, {'clang', 'gnu'}
    build_flag : str
    '''
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

    cxx_overwrite = ''
    if IS_OSX:
        if compiler == 'clang' or 'clang' in os.getenv('CXX', ''):
            cxx_overwrite = 'CXX="clang++ -stdlib=libstdc++"'
    print('cxx_overwrite flag', cxx_overwrite)

    cm = 'bash configure {build_flag} {compiler} {cxx_overwrite}'.format(
            build_flag=build_flag, compiler=compiler, cxx_overwrite=cxx_overwrite)

    print('build command: ', cm)
    # do not use subprocess to avoid split cxx_overwrite command
    os.system(cm)

    if IS_OSX:
        add_cpptraj_cxx_to_config('config.h', CPPTRAJ_CXX)

    subprocess.check_call('make libcpptraj -j4'.split())
    os.chdir(cwd)

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
    compiler, build_flag = get_compiler_and_build_flag()
    install_libcpptraj(compiler, build_flag)
