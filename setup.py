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
from scripts.base_setup import (add_openmp_flag, try_updating_libcpptraj)
from scripts.base_setup import CleanCommand, ISRELEASED, message_pip_need_cpptraj_home
from scripts.install_libcpptraj import DEFAULT_MAC_CCOMPILER, DEFAULT_MAC_CXXCOMPILER # clang

# python version >= 2.6
if sys.version_info < (2, 6):
    print('You must have at least Python 2.6 for pytraj\n')
    sys.exit(1)

amber_release = check_flag('--amber_release')
disable_openmp = check_flag('--disable-openmp')
openmp_flag = '-openmp' if not disable_openmp else ''
debug = check_flag('--debug')
use_phenix_python = check_flag('--phenix')
create_tar_file_for_release = True if 'sdist' in sys.argv else False
rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"
cpptraj_home = os.environ.get('CPPTRAJHOME', '')
use_pip = any('pip' in arg for arg in sys.argv)

if not cpptraj_home and use_pip:
    # if pip, require to set CPPTRAJHOME
    raise EnvironmentError(message_pip_need_cpptraj_home)

has_cpptraj_in_current_folder = os.path.exists("./cpptraj/")
phenix_python_lib = os.path.join(os.environ.get('PHENIX'),
                                 'base/lib/') if use_phenix_python else ''
pytraj_dir = os.path.abspath(os.path.dirname(__file__))
do_install, do_build = do_what(pytraj_dir)
cpptraj_dir, cpptraj_include, cpptraj_libdir, pytraj_inside_amber = get_include_and_lib_dir(rootname, cpptraj_home,
        has_cpptraj_in_current_folder, do_install, do_build, pytraj_dir, openmp_flag)
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
    print("adding debug info")
else:
    define_macros = []

# get INSTALLTYPE type from amber
installtype = os.environ.get("INSTALLTYPE", "")

if installtype:
    sys.argv.remove(installtype)

def read(fname):
    # must be in this setup file
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if not pytraj_inside_amber:
    if sys.platform == 'darwin':
        os.environ['CXX'] = DEFAULT_MAC_CXXCOMPILER
        os.environ['CC'] = DEFAULT_MAC_CCOMPILER
        # See which c++ lib we need to link to... sigh.
        import distutils.sysconfig as sc
        osxver = tuple(int(x) for x in
                       sc.get_config_var('MACOSX_DEPLOYMENT_TARGET').split('.') if x)
        if osxver < (10, 9):
            import platform
            minorosxver = int(platform.mac_ver()[0].split('.')[1])
            if minorosxver > 8:
                # OS X 10.8 and earlier do not understand this flag.
                extra_compile_args.extend(['-stdlib=libstdc++',
                                           '-mmacosx-version-min=%d.%d' % osxver])
                extra_link_args.extend(['-stdlib=libstdc++',
                                        '-mmacosx-version-min=%d.%d' % osxver])
else:
   # should use CXX and CC from config.h
   amberhome = os.environ.get('AMBERHOME', '')
   if not amberhome:
       raise EnvironmentError('must set AMBERHOME')

   configfile = amberhome + '/config.h'
   if not os.path.exists(configfile):
       raise OSError("must have config.h file")

   if sys.platform.startswith('darwin'):
       CC = DEFAULT_MAC_CCOMPILER
       CXX = DEFAULT_MAC_CXXCOMPILER
   elif sys.platform.startswith('linux'):
       CC='gcc'
       CXX='g++'
   else:
       pass

   with open(configfile) as fh:
       lines = fh.readlines()
       for line in lines:
           if line.startswith('CC='):
               CC = line.split('=', 1)[-1]
               break

       for line in lines:
           if line.startswith('CXX='):
               CXX = line.split('=', 1)[-1]
               break

   os.environ['CXX'] = CXX
   os.environ['CC'] = CC
   print('using CC={}, CXX={}'.format(CC, CXX))


pyxfiles, pxdfiles = get_pyx_pxd()

if not create_tar_file_for_release:
    if not libcpptraj_files:
        libcpptraj_files = try_updating_libcpptraj(cpptraj_home,
                do_install, do_build, has_cpptraj_in_current_folder, openmp_flag)
    print('libcpptraj_files', libcpptraj_files)

    try:
        output_openmp_check = subprocess.check_output(['nm', libcpptraj_files[0]]).decode().split('\n')
    except IndexError:
        print("Warning:  It seems that there is no libcpptraj. Please install it")
        sys.exit(1)

    libcpptraj_has_openmp = ([line for line in output_openmp_check if 'omp_get_num_threads' in line.lower()]  != [])
    if libcpptraj_has_openmp and sys.platform == 'darwin':
        raise OSError("pytraj does not (yet) support openmp in osx. Please recompile libcpptraj without openmp")

    if sys.platform == 'darwin' or sys.platform.startswith("win"):
        sys.stdout.write('does not support openmp on osx/win - disable\n')
        disable_openmp = True

    extra_compile_args, extra_link_args = add_openmp_flag(disable_openmp,
        libcpptraj_has_openmp, extra_compile_args, extra_link_args)

    if sys.platform.startswith('linux'):
        # set rpath
        sys.stdout.write('set rpath to {}\n'.format(cpptraj_libdir))
        extra_link_args.append('-Wl,-rpath={}'.format(cpptraj_libdir))
        extra_compile_args.append('-Wl,-rpath={}'.format(cpptraj_libdir))

    check_cpptraj_version(cpptraj_include, (4, 2, 8))

    pyxfiles, pxdfiles = get_pyx_pxd()

    if not do_clean and not ISRELEASED:
        if sys.platform.startswith("win"):
            cythonize(
                [pfile + '.pyx' for pfile in pyxfiles],
                compiler_directives=cython_directives,
            )
        else:
            cythonize(
                [pfile + '.pyx' for pfile in pyxfiles],
                nthreads=int(os.environ.get('NUM_THREADS', 4)),
                compiler_directives=cython_directives,
            )

    library_dirs = [cpptraj_libdir, ] if not use_phenix_python else [cpptraj_libdir, phenix_python_lib]

    if sys.platform.startswith('darwin') and use_pip:
        # ship with libcpptraj.dylib in pytraj/lib/
        try:
            os.mkdir('pytraj/lib')
        except OSError:
            pass

        shutil.copy('{}/libcpptraj.dylib'.format(cpptraj_libdir), 'pytraj/lib')
        os.system('install_name_tool -id @rpath/libcpptraj.dylib pytraj/lib/libcpptraj.dylib')
        library_dirs = ['pytraj/lib',]

    ext_modules = []
    for ext_name in pyxfiles:
        if need_cython:
            ext = ".pyx"
        else:
            ext = ".cpp"
        pyxfile = ext_name + ext

        # replace "/" by "." get module
        if "/" in ext_name:
            ext_name = ext_name.replace("/", ".")

        sources = [pyxfile]
        extmod = Extension(ext_name,
                           sources=sources,
                           libraries=['cpptraj'],
                           language='c++',
                           library_dirs=library_dirs,
                           define_macros=define_macros,
                           include_dirs=[cpptraj_include, pytraj_home],
                           extra_compile_args=extra_compile_args,
                           extra_link_args=extra_link_args)
        ext_modules.append(extmod)

else:
    # just need to create tar file for sdist
    ext_modules = []

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
sample_data = ["datafiles/ala3/Ala3.*",
               "datafiles/tz2/tz2.*",
               "datafiles/rna.pdb",
               "datafiles/trpcage/trpcage*",
               "datafiles/remd_ala2/*",
               "datafiles/dpdp/DPDP*"]

datalist = pxdfiles + sample_data

if sys.platform.startswith('darwin') and use_pip:
    datalist.append('lib/libcpptraj.dylib')

def build_func(ext_modules):
    return setup(
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


if __name__ == "__main__":
    build_func(ext_modules)
