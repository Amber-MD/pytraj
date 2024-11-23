#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import subprocess
import shutil
sys.path.append('scripts')
from check_openmp import get_openmp_flag
from find_lib import find_lib

DEFAULT_MAC_BUILD = '-shared -macAccelerate -noarpack'
DEFAULT_MAC_CCOMPILER = 'clang'
DEFAULT_MAC_CXXCOMPILER = 'clang++'
CPPTRAJ_CXX = ' '
IS_OSX = sys.platform.startswith('darwin')

if IS_OSX:
    import distutils.sysconfig as sc
    compiler = sc.get_config_var('CC') or ''
    if compiler == 'gcc':
        DEFAULT_MAC_CCOMPILER = "/usr/bin/gcc"
        DEFAULT_MAC_CXXCOMPILER = "/usr/bin/g++"

def add_cpptraj_cxx_to_config(fn, cpptraj_cxx=CPPTRAJ_CXX):
    with open(fn, 'r') as fconfig:
        lines = fconfig.readlines()
    with open('tmp.h', 'w') as fh:
        for idx, line in enumerate(lines):
            if line.startswith('CXX='):
                lines[idx] = ' '.join((line.strip(), cpptraj_cxx, '\n'))
        fh.write(''.join(lines))
    subprocess.check_call('mv tmp.h {}'.format(fn), shell=True)

def is_clang(cc):
    out = subprocess.check_output([cc, '--version']).decode()
    return 'clang' in out

def ensure_gnu():
    cc = os.getenv('CC', '')
    gcc_exe = cc if cc else 'gcc'
    if is_clang(gcc_exe):
        print('{} --version'.format(gcc_exe))
        print('{} here is actually clang compiler. Please export correct PATH for the real g++\n'.format(gcc_exe))
        print('Or export CXX and CC environments')
        print('e.g: ')
        print('    export CC=/usr/local/bin/gcc-5')
        print('    export CXX=/usr/local/bin/g++-5')
        sys.exit(1)

def get_compiler_and_build_flag():
    args = parse_args()
    openmp_flag = '-openmp' if args.openmp else ''
    if openmp_flag and not IS_OSX:
        assert get_openmp_flag(), 'your system must support openmp'
    install_type = args.install_type
    debug = '-debug' if args.debug else ''
    has_numpy = check_numpy()
    compiler, build_flag_ = determine_compiler_and_build_flag(has_numpy)
    amberhome = os.environ.get('AMBERHOME', '')
    amberlib = '-amberlib' if amberhome and args.amberlib else ''
    build_flag = construct_build_flag(build_flag_, amberlib, openmp_flag, debug)
    if install_type == 'github':
        print('install libcpptraj from github')
        subprocess.check_call('git clone https://github.com/Amber-MD/cpptraj'.split())
    else:
        print('install libcpptraj from current ./cpptraj folder')
    return compiler, build_flag

def check_numpy():
    try:
        import numpy as np
        return True
    except ImportError:
        return False

def determine_compiler_and_build_flag(has_numpy):
    os_dependent_default_compiler = 'clang' if IS_OSX else 'gnu'
    compiler_env = os.environ.get('COMPILER', '')
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
    prefix = sys.prefix
    if has_numpy and find_lib('openblas'):
        build_flag_ = '--with-blas={prefix} -openblas -noarpack'.format(prefix=prefix)
    elif has_numpy:
        build_flag_ = handle_numpy_without_openblas()
    else:
        build_flag_ = '-noarpack'
    zip_stuff = ' --with-bzlib={prefix} --with-zlib={prefix} '.format(prefix=prefix)
    build_flag_ += zip_stuff
    return compiler, build_flag_

def handle_numpy_without_openblas():
    import numpy as np
    try:
        blas_prefix = np.__config__.blas_opt_info['library_dirs'][0].strip('lib')
        lapack_prefix = np.__config__.lapack_opt_info['library_dirs'][0].strip('lib')
        return '-noarpack --with-blas={blas_prefix} --with-lapack={lapack_prefix}'.format(blas_prefix=blas_prefix, lapack_prefix=lapack_prefix)
    except (KeyError, IndexError):
        return '-noarpack'

def construct_build_flag(build_flag_, amberlib, openmp_flag, debug):
    if IS_OSX:
        build_flag = ' '.join((DEFAULT_MAC_BUILD, amberlib))
    else:
        build_flag = '--requires-flink ' + ' '.join(('-shared', build_flag_, amberlib, openmp_flag))
    build_flag = ' '.join((build_flag, debug, '-nosanderlib'))
    return build_flag

def install_libcpptraj(compiler='', build_flag='', n_cpus=4):
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
    cm = 'bash configure {build_flag} {compiler}'.format(build_flag=build_flag, compiler=compiler)
    print('configure command: ', cm)
    subprocess.check_call(cm, shell=True)
    if IS_OSX:
        add_cpptraj_cxx_to_config('config.h', CPPTRAJ_CXX)
    subprocess.check_call('make libcpptraj -j{}'.format(n_cpus).split())
    os.chdir(cwd)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='install options')
    parser.add_argument('-openmp', action='store_true', help='use openmp')
    parser.add_argument('-amberlib', action='store_true', help='use use amberlib if $AMBERHOME is set')
    parser.add_argument('-debug', action='store_true', help='debug')
    parser.add_argument('-j', default=4, help='n_cpus')
    parser.add_argument('install_type', default='', nargs='?', help='install_type in amber')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    compiler, build_flag = get_compiler_and_build_flag()
    install_libcpptraj(compiler, build_flag, args.j)
