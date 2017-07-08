'''pytraj: A python package binding to cpptraj program
'''

# install rule
# - if pytraj is inside $AMBERHOME, use libcpptraj.so in $AMBERHOME/lib and header file in cpptraj/src folder
# - if pytraj is outside $AMBERHOME
#    - check CPPTRAJ_LIBDIR and CPPTRAJ_HEADERDIR: if found, use those to install
#    - if not CPPTRAJ_LIBDIR, CPPTRAJ_HEADERDIR: check CPPTRAJHOME and found libcpptraj.so and header files in
#    CPPTRAJHOME/{lib, src}
#    - if not CPPTRAJHOME, find cpptraj folder in current folder
#    - if can not find cpptraj folder, do git clone from github

import os
import sys
import shutil

print('sys.argv', sys.argv)
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
from scripts.base_setup import (check_flag, write_version_py, get_version_info,
                                get_cpptraj_info, check_compile_cython,
                                check_cython, get_package_data, get_pyx_pxd)
from scripts.base_setup import (compiler_env_info, setenv_cc_cxx,
                                get_ext_modules)
from scripts.base_setup import CleanCommand, is_released

# python version >= 2.6
if sys.version_info < (2, 6):
    print('You must have at least Python 2.6 for pytraj\n')
    sys.exit(1)

try:
    if check_flag('--cythonize'):
        from Cython.Build import cythonize
        pyxfiles, _ = get_pyx_pxd()
        cythonize(
            [pfile + '.pyx' for pfile in pyxfiles], )
        sys.exit(0)
except ImportError:
    pass

amber_release = check_flag('--amber_release')
disable_openmp = check_flag('--disable-openmp')
use_amberlib = not check_flag('--disable-amberlib')
use_prebuilt_cythonized_files = check_flag('--use-pre-cythonized')
openmp_flag = '-openmp' if not disable_openmp else ''
debug = check_flag('-debug')
tarfile = True if 'sdist' in sys.argv else False
rootname = os.getcwd()
pytraj_src = os.path.join(rootname, "pytraj")
cpptraj_home = os.environ.get('CPPTRAJHOME', '')

compiler_env_info('Original compiler envs')
# danger
use_pip = (any('egg_info' in arg
               for arg in sys.argv) or any('pip' in arg for arg in sys.argv) or
           any('--no-deps' in arg for arg in sys.argv))

print('use_pip = {}, cpptraj_home = {}'.format(use_pip, cpptraj_home))
if use_pip and not cpptraj_home:
    if not os.path.exists('./cpptraj'):
        print(
            'Detect using pip, must set CPPTRAJHOME if there is no cpptraj in current folder'
        )
        sys.exit(1)

install_type = os.environ.get("INSTALLTYPE", "")
check_flag(install_type)
if install_type:
    print('install_type', install_type)

cpptraj_included = os.path.exists("./cpptraj/")
pytraj_home = os.path.abspath(os.path.dirname(__file__))

compile_c_extension = check_compile_cython(pytraj_home, use_pip=use_pip)
print('compile_c_extension =', compile_c_extension)

cpptraj_info = get_cpptraj_info(
    rootname=rootname,
    cpptraj_home=cpptraj_home,
    cpptraj_included=cpptraj_included,
    compile_c_extension=compile_c_extension,
    pytraj_home=pytraj_home,
    openmp_flag=openmp_flag,
    use_amberlib=use_amberlib)

libcpptraj_files = glob(os.path.join(cpptraj_info.lib_dir, 'libcpptraj') + '*')

write_version_py()
FULLVERSION, GIT_REVISION = get_version_info()
print(FULLVERSION)

# python setup.py clean
cmdclass = {'clean': CleanCommand}
need_cython, cmdclass, cythonize = check_cython(
    is_released,
    cmdclass,
    min_version='0.21',
    use_prebuilt_cythonized_files=use_prebuilt_cythonized_files)

if sys.platform.startswith('win'):
    extra_compile_args = []
    extra_link_args = []
else:
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
        'binding': True
    })
    define_macros = [
        ('CYTHON_TRACE', 1),
    ]
    print("adding debug info", cython_directives)
else:
    define_macros = []

setenv_cc_cxx(cpptraj_info.ambertools_distro, extra_compile_args,
              extra_link_args)

if not compile_c_extension:
    ext_modules = []
    possible_libcpptraj_files = []
else:
    ext_modules = get_ext_modules(
        cpptraj_info=cpptraj_info,
        pytraj_src=pytraj_src,
        compile_c_extension=compile_c_extension,
        is_released=is_released,
        need_cython=need_cython,
        cpptraj_included=cpptraj_included,
        libcpptraj_files=libcpptraj_files,
        openmp_flag=openmp_flag,
        use_amberlib=use_amberlib,
        cython_directives=cython_directives,
        Extension=Extension,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        define_macros=define_macros,
        use_pip=use_pip,
        tarfile=tarfile,
        use_prebuilt_cythonized_files=use_prebuilt_cythonized_files)
    if sys.platform.startswith('win'):
        possible_libcpptraj_files = [
            'libcpptraj.dll.a',
        ]
        for fn in libcpptraj_files:
            print("Copying {} to pytraj".format(fn))
            shutil.copy(fn, 'pytraj/')
    else:
        possible_libcpptraj_files = []

setup_args = {}
packages = [
    'pytraj',
    'pytraj.utils',
    'pytraj.builder',
    'pytraj.analysis',
    'pytraj.analysis.c_action',
    'pytraj.analysis.c_analysis',
    'pytraj.datasets',
    'pytraj.externals',
    'pytraj.trajectory',
    'pytraj.trajectory.c_traj',
    'pytraj.topology',
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
    'pytraj.serialize',
    'pytraj.sandbox',
    'pytraj.builder',
    'pytraj.testing',
]

if __name__ == "__main__":
    setup(
        name="pytraj",
        version=FULLVERSION,
        author="Hai Nguyen",
        url="https://github.com/Amber-MD/pytraj",
        packages=packages,
        description=
        """Python API for cpptraj: a data analysis package for biomolecular simulation""",
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
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Cython',
            'Programming Language :: C',
            'Programming Language :: C++',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
        ],
        ext_modules=ext_modules,
        package_data={
            'pytraj': get_package_data(use_pip) + possible_libcpptraj_files
        },
        cmdclass=cmdclass, )
