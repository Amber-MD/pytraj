#!/usr/bin/env python
from __future__ import print_function
import os, sys
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
    import numpy
    has_numpy = True
except ImportError:
    has_numpy = False

if has_numpy and find_lib('openblas'):
   prefix = sys.base_prefix
   # likely having openblas?
   build_flag = '--with-netcdf={prefix} --with-blas={prefix} --with-bzlib={prefix} --with-zlib={prefix} -openblas -noarpack'.format(prefix=prefix)
else:
   # user gets lucky?
   build_flag = ''

print('build_flag', build_flag)

cwd = os.getcwd()

compiler = os.environ.get('COMPILER', 'gnu')
amberhome = os.environ.get('AMBERHOME', '')
amberlib = '-amberlib' if amberhome else ''

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
              amberlib=amberlib,
              openmp=openmp_flag,
              build_flag=build_flag)

os.system('bash configure -shared {build_flag} {openmp} {amberlib} {compiler} || exit 1'.format(**config))

os.system('make libcpptraj -j8 || exit 1')
os.chdir(cwd)

print("make sure to 'export CPPTRAJHOME=$CPPTRAJHOME'"
"and 'export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:\$LD_LIBRARY_PATH'"
"then 'python ./setup.py install'")
