#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import subprocess
from subprocess import CalledProcessError
sys.path.append('./scripts')
from check_openmp import get_openmp_flag
from find_lib import find_lib

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


cwd = os.getcwd()

compiler = os.environ.get('COMPILER', 'gnu')
amberhome = os.environ.get('AMBERHOME', '')
amberlib = '-amberlib' if amberhome else ''

if has_numpy and find_lib('openblas'):
    prefix = sys.prefix
    # likely having openblas?
    build_flag = '--with-netcdf={prefix} --with-blas={prefix} \
                  --with-bzlib={prefix} --with-zlib={prefix} \
                  -openblas -noarpack'.format(prefix=prefix)
elif has_numpy:
    try:
        blas_prefix = np.__config__.blas_opt_info['library_dirs'][0].strip('lib')
        lapack_prefix = np.__config__.lapack_opt_info['library_dirs'][0].strip('lib')
        build_flag = '-noarpack --with-blas={blas_prefix} --with-lapack={lapack_prefix}'.format(blas_prefix=blas_prefix, lapack_prefix=lapack_prefix)
    except (KeyError, IndexError):
        build_flag = '-noarpack'
else:
    # user gets lucky?
    build_flag = '-noarpack'

build_flag = ' '.join((build_flag, amberlib, openmp_flag))

if install_type == 'github':
    print('install libcpptraj from github')
    os.system('git clone https://github.com/Amber-MD/cpptraj')
else:
    print('install libcpptraj from current ./cpptraj folder')

try:
    os.chdir('./cpptraj')
except FileNotFoundError:
    raise FileNotFoundError('please try ./scripts/install_cpptraj.py github')
os.environ['CPPTRAJHOME'] = os.getcwd()

try:
    os.mkdir('./lib')
except FileExistsError:
    pass

# turn off openmp. need to install pytraj with openmp too. Too complicated.
config = dict(compiler=compiler,
              build_flag=build_flag)

try:
    # assume that user has all required softwares
    cmd = 'bash configure -shared {openmp_flag} {amberlib} -noarpack {compiler}'.format(openmp_flag=openmp_flag, compiler=compiler, amberlib=amberlib)
    print('cmd', cmd)
    subprocess.check_call(cmd, shell=True)
except CalledProcessError:
    print('build_flag = ', build_flag)
    os.system('bash configure -shared {build_flag} {compiler} || exit 1'.format(**config))

os.system('make libcpptraj -j8 || exit 1')
os.chdir(cwd)

print("make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
      "and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
      "then 'python ./setup.py install'")
