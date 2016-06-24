'''install rule
- if pytraj is inside $AMBERHOME, use libcpptraj.so in $AMBERHOME/lib and header file in cpptraj/src folder
- if pytraj is outside $AMBERHOME
    - check CPPTRAJ_LIBDIR and CPPTRAJ_HEADERDIR: if found, use those to install
    - if not CPPTRAJ_LIBDIR, CPPTRAJ_HEADERDIR: check CPPTRAJHOME and found libcpptraj.so and header files in
    CPPTRAJHOME/{lib, src}
    - if not CPPTRAJHOME, find cpptraj folder in current folder
    - if can not find cpptraj folder, do git clone from github
'''

import os
import sys
import subprocess
import shutil
try:
    # for amber
    sys.argv.remove('--no-setuptools')
    from distutils.core import setup
    from distutils.extension import Extension
except ValueError:
    try:
        from setuptools import setup, Extension
    except ImportError:
        from distutils.core import setup
        from distutils.extension import Extension
from glob import glob


# local import
from scripts.base_setup import (check_flag, check_cpptraj_version, write_version_py, get_version_info,
                                get_pyx_pxd, get_include_and_lib_dir, do_what, check_cython)
from scripts.base_setup import (add_openmp_flag, try_updating_libcpptraj, setenv_cc_cxx, get_ext_modules)
from scripts.base_setup import CleanCommand, ISRELEASED, message_pip_need_cpptraj_home
from scripts.install_libcpptraj import DEFAULT_MAC_CCOMPILER, DEFAULT_MAC_CXXCOMPILER # clang

# python version >= 2.6
if sys.version_info < (2, 6):
    print('You must have at least Python 2.6 for pytraj\n')
    sys.exit(1)

amber_release = check_flag('--amber_release')
disable_openmp = check_flag('--disable-openmp')
use_amberlib = not check_flag('--disable-amberlib')
openmp_flag = '-openmp' if not disable_openmp else ''
debug = check_flag('-debug')
tarfile = True if 'sdist' in sys.argv else False
rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"
cpptraj_home = os.environ.get('CPPTRAJHOME', '')
use_pip = any('pip' in arg for arg in sys.argv)

if not cpptraj_home and use_pip:
    # if pip, require to set CPPTRAJHOME
    raise EnvironmentError(message_pip_need_cpptraj_home)

cpptraj_included = os.path.exists("./cpptraj/")
pytraj_dir = os.path.abspath(os.path.dirname(__file__))
do_install, do_build = do_what(pytraj_dir)
cpptraj_dir, cpptraj_include_dir, cpptraj_libdir, ambertools_distro = get_include_and_lib_dir(rootname, cpptraj_home,
        cpptraj_included, do_install, do_build, pytraj_dir, openmp_flag)
libcpptraj_files = glob(os.path.join(cpptraj_libdir, 'libcpptraj') + '*')
do_clean = (len(sys.argv) == 2 and 'clean' in sys.argv)

write_version_py()
FULLVERSION, GIT_REVISION = get_version_info()
print(FULLVERSION)

# python setup.py clean
cmdclass = {'clean': CleanCommand}
need_cython, cmdclass, cythonize  = check_cython(ISRELEASED, cmdclass, min_version='0.21')

extra_compile_args = ['-O0', '-ggdb']
extra_link_args = ['-O0', '-ggdb']

cython_directives = {
    'embedsignature': True,
    'boundscheck': False,
    'wraparound': False,
}

if debug:
    cython_directives.update({
        'profile': True,
        'linetrace': True,
        'binding': True})
    define_macros = [('CYTHON_TRACE', 1), ]
    print("adding debug info", cython_directives)
else:
    define_macros = []

# get INSTALLTYPE type from amber
installtype = os.environ.get("INSTALLTYPE", "")

if installtype:
    print('install_type', install_type)
    sys.argv.remove(installtype)

setenv_cc_cxx(ambertools_distro)

ext_modules = get_ext_modules(cpptraj_home,
                cpptraj_libdir,
                cpptraj_include_dir,
                pytraj_home,
                do_install,
                do_build,
                do_clean,
                ISRELEASED,
                cpptraj_included,
                libcpptraj_files,
                openmp_flag,
                use_amberlib,
                extra_compile_args=[],
                extra_link_args=[],
                define_macros=[],
                tarfile=False)

setup_args = {}
packages = [
    'pytraj',
    'pytraj.utils',
    'pytraj.c_action',
    'pytraj.c_analysis',
    'pytraj.datasets',
    'pytraj.externals',
    'pytraj.c_traj',
    'pytraj.datafiles',
    'pytraj.datafiles.ala3',
    'pytraj.datafiles.tz2',
    'pytraj.datafiles.dpdp',
    'pytraj.datafiles.trpcage',
    'pytraj.datafiles.remd_ala2',
    'pytraj.plot',
    'pytraj.math',
    'pytraj.core',
    'pytraj.parallel',
    'pytraj.cluster',
    'pytraj.visualization',
    'pytraj.sandbox',
]

pylen = len('pytraj') + 1
sample_datafiles  = ["datafiles/ala3/Ala3.*",
               "datafiles/tz2/tz2.*",
               "datafiles/rna.pdb",
               "datafiles/trpcage/trpcage*",
               "datafiles/remd_ala2/*",
               "datafiles/dpdp/DPDP*"]

jsfiles = ['utils/progress-circle/css/*css',
      'utils/progress-circle/*js',]

datalist = pxdfiles + sample_datafiles + jsfiles

if sys.platform.startswith('darwin') and use_pip:
    datalist.append('lib/libcpptraj.dylib')

if __name__ == "__main__":
    setup(
        name="pytraj",
        version=FULLVERSION,
        author="Hai Nguyen",
        url="https://github.com/Amber-MD/pytraj",
        packages=packages,
        description="""Python API for cpptraj: a data analysis package for biomolecular simulation""",
        license="GPL v3",
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Operating System :: Unix',
            'Operating System :: MacOS',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Cython',
            'Programming Language :: C',
            'Programming Language :: C++',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
        ],
        ext_modules=ext_modules,
        package_data={'pytraj': datalist},
        cmdclass=cmdclass,)
