#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import subprocess
import shutil
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


def is_clang(cc):
    # both cpptraj and pytraj give priority for CXX and CC environments
    # check them first.
    out = subprocess.check_output([cc, '--version']).decode()
    return 'clang' in out


def ensure_gnu():
    # both cpptraj and pytraj give priority for CXX and CC environments
    # check them first.
    cc = os.getenv('CC', '')
    gcc_exe = cc if cc else 'gcc'

    if is_clang(gcc_exe):
        print('{} --version'.format(gcc_exe))
        print(
            '{} here is actually clang compiler. Please export correct PATH for the real g++\n'.
            format(gcc_exe))
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
    compiler_env = os.environ.get('COMPILER',
                                  '')  # intel | pgi | clang | cray?
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
        print('has_numpy but can not find openblas')
        print('try blas and lapack')
        try:
            blas_prefix = np.__config__.blas_opt_info['library_dirs'][0].strip(
                'lib')
            lapack_prefix = np.__config__.lapack_opt_info['library_dirs'][
                0].strip('lib')
            build_flag_ = (
                '-noarpack --with-blas={blas_prefix} --with-lapack={lapack_prefix}'
                .format(blas_prefix=blas_prefix, lapack_prefix=lapack_prefix))
        except (KeyError, IndexError):
            build_flag_ = '-noarpack'
    else:
        # user gets lucky?
        build_flag_ = '-noarpack'

    zip_stuff = ' --with-bzlib={prefix} --with-zlib={prefix} '.format(
        prefix=prefix)
    build_flag_ += zip_stuff

    if IS_OSX:
        build_flag = ' '.join((DEFAULT_MAC_BUILD, amberlib))
    else:
        build_flag = '--requires-flink ' + ' '.join(('-shared', build_flag_, amberlib, openmp_flag))
        # FIXME: remove --requires-flink if https://github.com/Amber-MD/cpptraj/issues/585 is resolved

    build_flag = ' '.join((build_flag, debug, '-nosanderlib'))

    if install_type == 'github':
        print('install libcpptraj from github')
        subprocess.check_call(
            'git clone https://github.com/Amber-MD/cpptraj'.split())
    else:
        print('install libcpptraj from current ./cpptraj folder')
    return compiler, build_flag


def install_libcpptraj(compiler='', build_flag='', n_cpus=4):
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

    if sys.platform.startswith('win'):
        _install_libcpptraj_win_msys2()
    else:
        cxx_overwrite = ''

        # We don't need this anymore?
        # if IS_OSX:
        #     if ((compiler == 'clang' or 'clang' in os.getenv('CXX', '')) or
        #         (os.getenv('CXX') and is_clang(os.getenv('CXX')))):
        #         # cxx_overwrite = 'CXX="clang++ -stdlib=libstdc++"'
        print('cxx_overwrite flag', cxx_overwrite)

        cm = 'bash configure {build_flag} {compiler} {cxx_overwrite}'.format(
            build_flag=build_flag,
            compiler=compiler,
            cxx_overwrite=cxx_overwrite)

        print('configure command: ', cm)
        # do not use subprocess to avoid split cxx_overwrite command
        os.system(cm)

        if IS_OSX:
            add_cpptraj_cxx_to_config('config.h', CPPTRAJ_CXX)

        subprocess.check_call('make libcpptraj -j{}'.format(n_cpus).split())
        os.chdir(cwd)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='install options')
    parser.add_argument('-openmp', action='store_true', help='use openmp')
    parser.add_argument(
        '-amberlib',
        action='store_true',
        help='use use amberlib if $AMBERHOME is set')
    parser.add_argument('-debug', action='store_true', help='debug')
    parser.add_argument('-j', default=4, help='n_cpus')
    parser.add_argument(
        'install_type', default='', nargs='?', help='install_type in amber')
    args = parser.parse_args()
    return args


def _install_libcpptraj_win_msys2():
    PREFIX = '/usr/local/'
    PREFIX2 = '/mingw64/'
    # assume you do all the setup
    command = """
    sh configure --with-netcdf={PREFIX} \
                 --with-blas={PREFIX2} \
                 --with-arpack={PREFIX2} \
                 --with-readline={PREFIX2} \
                 -openblas \
                 -shared \
                 -windows \
                 gnu
    """.format(
        PREFIX=PREFIX, PREFIX2=PREFIX2).strip()
    subprocess.check_call(command, shell=True)
    subprocess.check_call('make libcpptraj -j2', shell=True)
    # will create libcpptraj.dll.a
    # shutil.copy('lib/libcpptraj.so', 'lib/cpptraj.lib')


if __name__ == '__main__':
    if sys.platform.startswith('win'):
        install_libcpptraj()
    else:
        args = parse_args()
        compiler, build_flag = get_compiler_and_build_flag()
        install_libcpptraj(compiler, build_flag, args.j)
