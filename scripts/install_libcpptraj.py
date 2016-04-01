#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import subprocess
from subprocess import CalledProcessError
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
    os.system('mv tmp.h {}'.format(fn))


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

    # better name for CPPTRAJ_COMPILER_OPTION?
    # cpptraj: ./configure gnu
    # e.g: CPPTRAJ_COMPILER_OPTION=gnu python ./scripts/install_libcpptraj.py
    cpptraj_compiler_option = os.environ.get('CPPTRAJ_COMPILER_OPTION', 'gnu')  # intel | pgi | clang | cray?
    amberhome = os.environ.get('AMBERHOME', '')
    amberlib = '-amberlib' if amberhome else ''

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
        build_flag = DEFAULT_MAC_BUILD
    else:
        build_flag = ' '.join(('-shared', build_flag_, amberlib, openmp_flag))

    if install_type == 'github':
        print('install libcpptraj from github')
        os.system('git clone https://github.com/Amber-MD/cpptraj')
    else:
        print('install libcpptraj from current ./cpptraj folder')
    return cpptraj_compiler_option, build_flag


def fix_rpath():
    if IS_OSX:
        os.system('(cd lib && install_name_tool -id `pwd`/libcpptraj.dylib libcpptraj.dylib)')
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

    cm = 'bash configure {build_flag} {compiler} {cxx_overwrite}|| exit 1'.format(
            build_flag=build_flag, compiler=cpptraj_compiler_option, cxx_overwrite=cxx_overwrite)

    print('build command: ', cm)
    os.system(cm)

    if IS_OSX:
        add_CPPTRAJ_CXX_to_config('config.h', CPPTRAJ_CXX)

    os.system('make libcpptraj -j8 || exit 1')
    fix_rpath()
    os.chdir(cwd)

    print("make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
          "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
          "then 'python ./setup.py install'")

if __name__ == '__main__':
    cpptraj_compiler_option, build_flag = get_compiler_and_build_flag()
    install_libcpptraj(cpptraj_compiler_option, build_flag)

