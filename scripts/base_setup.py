#!/usr/bin/env python

import os
import sys
import time
import shutil
import subprocess
from distutils.command.clean import clean as Clean
from subprocess import CalledProcessError
from glob import glob

if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins

MAJOR = 1
MINOR = 0
MICRO = 2
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


message_cython = """
Warning:  pytraj was not installed !

Building from source requires cython >= 0.21

Either try:
    conda install cython
    pip install cython --upgrade

We suggest to install Anaconda suite, which has more than 300 python packages (including
the most updated cython)

    http://conda.pydata.org/docs/download.html)

"""
message_pip_need_cpptraj_home = """

installing from pip require to pre-install libcpptraj and to set CPPTRAJHOME

An example of installing libcpptraj:

$ git clone https://github.com/Amber-MD/cpptraj/
$ cd cpptraj
$ export CPPTRAJHOME=`pwd`
$ ./configure -shared -openmp gnu
$ make libcpptraj -j8

"""

message_auto_install = """
Can not find cpptraj header and libcpptraj files.
We're trying to download and build libcpptraj for you, would take about 5-10 minutes.
You can check ./cpptraj/ folder after installation.

To avoid auto-installation
--------------------------
You should set CPPTRAJHOME or installing ./cpptraj/ in this current folder.

If you want to manually install `libcpptraj`, you can download cpptraj
development version from here: https://github.com/Amber-MD/cpptraj

$ git clone https://github.com/Amber-MD/cpptraj/
$ cd cpptraj
$ export CPPTRAJHOME=`pwd`
$ ./configure -shared -openmp gnu
$ make libcpptraj

and then go back to pytraj folder:
python setup.py install
"""

message_openmp_cpptraj = '''
libcpptraj was detected to be installed with openmp.
You can not use --disable-openmp flag with pytraj
'''

message_serial_cpptraj = '''
libcpptraj was NOT installed with openmp. You can recompile it with -openmp flag or
disable openmp install in pytraj by adding --disable-openmp

Note: For osx users, pytraj uses clang for compiling cython extension and '--disable-openmp' flag
must be specified. If experienced users want to hack, please check setup.py file.

Example:
    - Turn off openmp in pytraj: python setup.py install --disable-openmp
    - Turn ON openmp in cpptraj: ./configure -shared -openmp gnu && make libcpptraj -j4
    - Turn ON openmp in cpptraj and using netcdf, lapack, ... in $AMBERHOME (if having
      one):
          ./configure -shared -openmp -amberlib gnu && make libcpptraj -j4
'''

message_after_sucessful_install = '''
make sure to add {0} to your LD_LIBRARY_PATH
    - Example: export LD_LIBRARY_PATH={1}:$LD_LIBRARY_PATH

    - Notes: you can move `libcpptraj.so` to any where you want, just properly add it to $LD_LIBRARY_PATH

Run test:
    - simple (few seconds): python ./run_tests.py simple
    - full (5-10 minutes): python run_tests.py
'''

def check_flag(key):
    try:
        sys.argv.remove(key)
        return True
    except ValueError:
        return False

def check_cpptraj_version(header_dir, version=(4, 2, 9)):
    vfile = os.path.join(header_dir, 'Version.h')
    with open(vfile) as fh:
        for line in fh.readlines():
            if line.startswith('#define CPPTRAJ_INTERNAL_VERSION'):
                break
        int_version = tuple(int(i) for i in line.split()[-1].strip('"').replace('V', '').split('.'))
        if int_version < version:
            sys.stderr.write('must have cpptraj version >= {}\n'.format(version))
            sys.exit(0)


def remind_export_LD_LIBRARY_PATH(build_tag, libdir, pytraj_inside_amber):
    if build_tag:
        if not pytraj_inside_amber:
            from scripts.acsii_art import batman
            print("")
            print("")
            print(batman)
            libdir = os.path.abspath(libdir)
            print(message_after_sucessful_install.format(libdir, libdir))
            print("")
        else:
            # pytraj is a part of Amber
            print('make sure to `source $AMBERHOME/amber.sh` (if using bash) '
                  'or `source $AMBERHOME/amber.csh` if using csh')
    else:
        print("not able to install pytraj")

# Return the git revision as a string
# git_version, get_version_info, write_version_py  was lightly adapted from numpy package
# http://www.numpy.org/
# Numpy is licensed under the BSD license,


def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

# This is a bit hackish: we are setting a global variable so that the main
# numpy __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.  While ugly, it's
# a lot more robust than what was previously being used.
builtins.__PYTRAJ_SETUP__ = True


def get_version_info():
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of numpy.version messes up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('pytraj/__version__.py'):
        # must be a source distribution, use existing version file
        try:
            from pytraj.__version__ import git_revision as GIT_REVISION
        except ImportError:
            # FIXME pytraj
            # raise ImportError("Unable to import git_revision. Try removing " \
            #                  "numpy/version.py and the build directory " \
            #                  "before building.")
            pass
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '-beta0'

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='pytraj/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM PYTRAJ SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

def do_what(PYTRAJ_DIR):
    # this checking should be here, after checking openmp and other stuff
    if len(sys.argv) == 2 and sys.argv[1] == 'install':
        do_install = True
    elif len(sys.argv) == 3 and sys.argv[1] == 'install' and os.path.join('AmberTools',
                                                                          'src') in PYTRAJ_DIR:
        # install pytraj in $AMBERHOME
        # do not use pytraj_inside_amber here in we call `do_what()` before calling get_include_and_lib_dir()
        # don't mess this up
        # $(PYTHON) setup.py install $(PYTHON_INSTALL)
        do_install = True
    else:
        do_install = False

    if len(sys.argv) == 2 and sys.argv[1] == 'build':
        do_build = True
    else:
        do_build = False
    return do_install, do_build

def install_libcpptraj(openmp_flag, from_github=False):
    github = 'github' if from_github else ''
    options = dict(github=github,
                   openmp_flag=openmp_flag)
    cmd = "python scripts/install_libcpptraj.py {github} {openmp_flag}".format(**options)
    print('command = ', cmd)
    subprocess.check_call(cmd, shell=True)

def try_updating_libcpptraj(cpptraj_home,
                            do_install,
                            do_build, 
                            has_cpptraj_in_current_folder,
                            openmp_flag):
    if cpptraj_home:
        raise ValueError(
            '$CPPTRAJHOME exists but there is no libcpptraj in $CPPTRAJHOME/lib \n'
            'There are two solutions: \n'
            '1. unset CPPTRAJHOME and `python setup.py install` again. We will install libcpptraj for you. \n'
            '2. Or you need to install libcpptraj in $CPPTRAJHOME/lib \n')
    else:
        if do_install or do_build:
            if has_cpptraj_in_current_folder:
                print(
                    'can not find libcpptraj but found ./cpptraj folder, trying to reinstall it to ./cpptraj/lib/ \n')
                time.sleep(3)
                try:
                    cpptraj_dir = './cpptraj/'
                    cpptraj_libdir = cpptraj_dir + '/lib/'
                    install_libcpptraj(openmp_flag, from_github=False)
                    return glob(os.path.join(cpptraj_libdir, 'libcpptraj') + '*')
                except CalledProcessError:
                    print(
                        'can not install libcpptraj, you need to install it manually \n')
                    sys.exit(0)
            else:
                print('can not find libcpptraj in $CPPTRAJHOME/lib. '
                                 'You need to install ``libcpptraj`` manually. ')
                sys.exit(0)


def add_openmp_flag(disable_openmp,
                    libcpptraj_has_openmp,
                    extra_compile_args,
                    extra_link_args):
    if disable_openmp:
        if libcpptraj_has_openmp:
            raise ValueError(message_openmp_cpptraj)
        else:
            return (extra_compile_args, extra_link_args)

    else:
        if not libcpptraj_has_openmp:
            raise ValueError(message_serial_cpptraj)
        # make copy
        return (extra_compile_args[:] + ["-fopenmp",], extra_link_args[:] + ["-fopenmp",])

def get_pyx_pxd():
    pxd_include_dirs = [
        directory for directory, dirs, files in os.walk('pytraj') if '__init__.pyx'
        in files or '__init__.pxd' in files or '__init__.py' in files
    ]
    
    pxd_include_patterns = [p + '/*.pxd' for p in pxd_include_dirs]
    
    pyxfiles = []
    for p in pxd_include_dirs:
        pyxfiles.extend([ext.split(".")[0] for ext in glob(p + '/*.pyx')
                         if '.pyx' in ext])
    pxdfiles = [p.replace("pytraj/", "") for p in pxd_include_patterns]
    return pyxfiles, pxdfiles

def check_cython(ISRELEASED, cmdclass, min_version='0.21'):
    if ISRELEASED:
        # ./devtools/mkrelease
        need_cython = False
        cythonize = None
    else:
        try:
            import Cython
            from Cython.Distutils import build_ext
            from Cython.Build import cythonize
            need_cython = True
            cmdclass['build_ext'] = build_ext
            if Cython.__version__ < min_version:
                print(message_cython)
                sys.exit(0)
        except ImportError:
            print(message_cython)
            sys.exit(0)
    return need_cython, cmdclass, cythonize


def get_include_and_lib_dir(rootname, cpptrajhome, has_cpptraj_in_current_folder, do_install, do_build, PYTRAJ_DIR, openmp_flag):
    # check if has environment variables
    CPPTRAJ_LIBDIR = os.environ.get('CPPTRAJ_LIBDIR', '')
    CPPTRAJ_HEADERDIR = os.environ.get('CPPTRAJ_HEADERDIR', '')
    if os.path.join('AmberTools', 'src') in PYTRAJ_DIR:
        # install pytraj inside AMBER
        AMBERHOME = os.environ.get('AMBERHOME', '')
        if not AMBERHOME:
            raise EnvironmentError('must set AMBERHOME if you want to install pytraj '
                                   'inside AMBER')
        # overwrite CPPTRAJ_HEADERDIR, CPPTRAJ_LIBDIR
        CPPTRAJ_LIBDIR = os.path.join(AMBERHOME, 'lib')
        CPPTRAJ_HEADERDIR = os.path.join(AMBERHOME, 'AmberTools', 'src', 'cpptraj', 'src')

        pytraj_inside_amber = True
    else:
        pytraj_inside_amber = False

    if CPPTRAJ_LIBDIR and CPPTRAJ_HEADERDIR:
        cpptraj_include = CPPTRAJ_HEADERDIR
        libdir = CPPTRAJ_LIBDIR
        cpptraj_dir = ''
    else:
        if cpptrajhome:
            # use libcpptraj and header files in CPPTRAJHOME (/lib, /src)
            cpptraj_dir = cpptrajhome
            cpptraj_include = cpptraj_dir + "/src/"
            libdir = cpptrajhome + "/lib/"
        elif has_cpptraj_in_current_folder:
            cpptraj_dir = os.path.abspath("./cpptraj/")
            cpptraj_include = cpptraj_dir + "/src/"
            libdir = cpptraj_dir + "/lib/"
        else:

            if do_install or do_build:
                print(message_auto_install)
                for i in range(0, 3):
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    time.sleep(1)
                try:
                    install_libcpptraj(openmp_flag, from_github=True)
                except CalledProcessError:
                    print(
                        'can not install libcpptraj, you need to install it manually \n')
                    sys.exit(0)
            cpptraj_dir = os.path.join(rootname, "cpptraj")
            cpptraj_include = os.path.join(cpptraj_dir, 'src')
            libdir = os.path.join(cpptraj_dir, 'lib')
    return cpptraj_dir, cpptraj_include, libdir, pytraj_inside_amber


# CleanCommand was copied and lightly adapted from scikit-learn package
# https://github.com/scikit-learn/scikit-learn
# New BSD License

# Custom clean command to remove build artifacts
class CleanCommand(Clean):
    description = "Remove build artifacts from the source tree"

    def run(self):
        Clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        for dirpath, dirnames, filenames in os.walk('pytraj'):
            for filename in filenames:
                if (filename.endswith('.so') or filename.endswith('.pyd')
                        or filename.endswith('.dll')
                        or filename.endswith('.pyc')):
                    os.unlink(os.path.join(dirpath, filename))
            for dirname in dirnames:
                if dirname == '__pycache__':
                    shutil.rmtree(os.path.join(dirpath, dirname))
