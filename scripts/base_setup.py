#!/usr/bin/env python

import os
import sys
import time
import shutil
import subprocess
from glob import glob
from subprocess import CalledProcessError
from distutils.command.clean import clean as Clean

if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins

MAJOR = 1
MINOR = 1
MICRO = 0
is_released = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

DEFAULT_MAC_CCOMPILER = 'clang'
DEFAULT_MAC_CXXCOMPILER = 'clang++'


message_cython = '''
Warning:  pytraj was not installed !

Building from source requires cython >= 0.21

Either try:
    conda install cython
    pip install cython --upgrade

We suggest to install Anaconda suite, which has more than 300 python packages (including
the most updated cython)

    http://conda.pydata.org/docs/download.html)

'''

message_pip_need_cpptraj_home = '''

installing from pip for python 2.7 in OSX
requires to pre-install libcpptraj and to set CPPTRAJHOME

An example of installing libcpptraj:

$ git clone https://github.com/Amber-MD/cpptraj/
$ cd cpptraj
$ export CPPTRAJHOME=`pwd`
$ ./configure -shared -openmp gnu
$ make libcpptraj -j8

How to install pytraj within seconds without installing libcpptraj by yourself?

Either:
    - Use python >= 3.4 (OSX): pip install pytraj
    - Use conda (Linux, osx): conda install -c ambermd pytraj
'''

message_auto_install = '''
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
'''

message_openmp_cpptraj = '''
libcpptraj was detected to be installed with openmp.
You can not use --disable-openmp flag with pytraj
'''

message_serial_cpptraj = '''
libcpptraj was NOT installed with openmp. You can recompile it with -openmp flag or
disable openmp install in pytraj by adding --disable-openmp

Note: For osx users, pytraj uses clang for compiling cython extension and '--disable-openmp' flag
must be specified. If experienced users want to hack, please check setup.py file.

Examples:
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
            sys.exit(1)


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

    if not is_released:
        FULLVERSION += '.dev0'

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='pytraj/version.py'):
    cnt = '''
# THIS FILE IS GENERATED FROM PYTRAJ SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
'''
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(is_released)})
    finally:
        a.close()

def check_compile_cython(pytraj_home, use_pip=False):
    compile_c_extension = False
    # this checking should be here, after checking openmp and other stuff
    if '--help' in sys.argv or '-h' in sys.argv or '--help-commands' in sys.argv:
        compile_c_extension = False
    elif use_pip:
        compile_c_extension = True
    else:
        if ('install' in sys.argv or
            'build' in sys.argv or
            'build_ext' in sys.argv):
            compile_c_extension = True
    return compile_c_extension

def install_libcpptraj(openmp_flag, from_github=False, use_amberlib=True):
    '''If AMBERHOME is set and amberlib is True, use -amberlib for libcpptraj
    '''
    github = 'github' if from_github else ''
    amberlib_flag = '-amberlib' if use_amberlib else ''
    options = dict(github=github,
                   openmp_flag=openmp_flag,
                   amberlib_flag=amberlib_flag)
    cmd = "python scripts/install_libcpptraj.py {github} {openmp_flag} {amberlib_flag}".format(**options)
    print('command = ', cmd)
    subprocess.check_call(cmd, shell=True)

def try_updating_libcpptraj(cpptraj_home,
                            compile_c_extension,
                            cpptraj_included,
                            openmp_flag,
                            use_amberlib):
    if cpptraj_home:
        raise ValueError(
            '$CPPTRAJHOME exists but there is no libcpptraj in $CPPTRAJHOME/lib \n'
            'There are two solutions: \n'
            '1. unset CPPTRAJHOME and `python setup.py install` again. We will install libcpptraj for you. \n'
            '2. Or you need to install libcpptraj in $CPPTRAJHOME/lib \n')
    else:
        if compile_c_extension:
            print('compile_c_extension')
            if cpptraj_included:
                print(
                    'can not find libcpptraj but found ./cpptraj folder, trying to reinstall it to ./cpptraj/lib/ \n')
                time.sleep(3)
                try:
                    cpptraj_dir = './cpptraj/'
                    cpptraj_lib_dir = cpptraj_dir + '/lib/'
                    install_libcpptraj(openmp_flag, from_github=False, use_amberlib=use_amberlib)
                    return glob(os.path.join(cpptraj_lib_dir, 'libcpptraj') + '*')
                except CalledProcessError:
                    print(
                        'can not install libcpptraj, you need to install it manually \n')
                    sys.exit(1)
            else:
                print('can not find libcpptraj in $CPPTRAJHOME/lib. '
                                 'You need to install ``libcpptraj`` manually. ')
                sys.exit(1)


def add_openmp_flag(openmp_flag,
                    libcpptraj_has_openmp,
                    extra_compile_args,
                    extra_link_args):
    if not openmp_flag:
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

def check_cython(is_released, cmdclass, min_version='0.21'):
    if is_released:
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
                sys.exit(1)
        except ImportError:
            print(message_cython)
            sys.exit(1)
    return need_cython, cmdclass, cythonize


class CpptrajInfo(object):
    home_env = None
    dir = None
    lib_dir = None
    int_version = None
    ambertools_distro = False

def get_cpptraj_info(rootname,
                     cpptraj_home,
                     cpptraj_included,
                     compile_c_extension,
                     pytraj_home,
                     openmp_flag,
                     use_amberlib):

    cpptraj_info = CpptrajInfo()

    # check if has environment variables
    CPPTRAJ_LIBDIR = os.environ.get('CPPTRAJ_LIBDIR', '')
    CPPTRAJ_HEADERDIR = os.environ.get('CPPTRAJ_HEADERDIR', '')
    if os.path.join('AmberTools', 'src') in pytraj_home:
        # install pytraj inside AMBER
        AMBERHOME = os.environ.get('AMBERHOME', '')
        if not AMBERHOME:
            raise EnvironmentError('must set AMBERHOME if you want to install pytraj '
                                   'inside AMBER')
        # overwrite CPPTRAJ_HEADERDIR, CPPTRAJ_LIBDIR
        CPPTRAJ_LIBDIR = os.path.join(AMBERHOME, 'lib')
        CPPTRAJ_HEADERDIR = os.path.join(AMBERHOME, 'AmberTools', 'src', 'cpptraj', 'src')

        cpptraj_info.ambertools_distro = True
    else:
        cpptraj_info.ambertools_distro = False

    if CPPTRAJ_LIBDIR and CPPTRAJ_HEADERDIR:
        cpptraj_info.include_dir = CPPTRAJ_HEADERDIR
        cpptraj_info.lib_dir = CPPTRAJ_LIBDIR
        cpptraj_info.dir = ''
    else:
        if cpptraj_home:
            # use libcpptraj and header files in CPPTRAJHOME (/lib, /src)
            cpptraj_info.dir = cpptraj_home
            cpptraj_info.include_dir = cpptraj_info.dir + "/src/"
            cpptraj_info.lib_dir = cpptraj_info.dir + "/lib/"
        elif cpptraj_included:
            cpptraj_info.dir = os.path.abspath("./cpptraj/")
            cpptraj_info.include_dir = cpptraj_info.dir + "/src/"
            cpptraj_info.lib_dir = cpptraj_info.dir + "/lib/"
        else:

            if compile_c_extension:
                print(message_auto_install)
                for i in range(0, 3):
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    time.sleep(1)
                try:
                    install_libcpptraj(openmp_flag, from_github=True, use_amberlib=use_amberlib)
                except CalledProcessError:
                    print(
                        'can not install libcpptraj, you need to install it manually \n')
                    sys.exit(1)
            cpptraj_info.dir = os.path.join(rootname, "cpptraj")
            cpptraj_info.include_dir = os.path.join(cpptraj_info.dir, 'src')
            cpptraj_info.lib_dir = os.path.join(cpptraj_info.dir, 'lib')
    return cpptraj_info

def setenv_cc_cxx(ambertools_distro,
        extra_compile_args,
        extra_link_args):
    '''force pytraj and cpptraj to use the sample compiler if pytraj
    is distribued by AmberTools.
    '''
    if not ambertools_distro:
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
        print('pytraj is inside AMBERHOME')
        # should use CXX and CC from config.h
        amber_home = os.environ.get('AMBERHOME', '')
        if not amber_home:
            raise EnvironmentError('must set AMBERHOME')

        configfile = amber_home + '/config.h'
        if not os.path.exists(configfile):
            raise OSError("must have config.h file")

        # make default compiler first
        if sys.platform.startswith('darwin'):
            CC = DEFAULT_MAC_CCOMPILER
            CXX = DEFAULT_MAC_CXXCOMPILER
        elif sys.platform.startswith('linux'):
            CC='gcc'
            CXX='g++'
        else:
            pass

        # then parse $AMBERHOME/config.h
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

def get_ext_modules(cpptraj_info,
                    pytraj_src,
                    compile_c_extension,
                    is_released,
                    need_cython,
                    cpptraj_included,
                    libcpptraj_files,
                    openmp_flag,
                    use_amberlib,
                    cython_directives,
                    Extension,
                    extra_compile_args=[],
                    extra_link_args=[],
                    define_macros=[],
                    use_pip=False,
                    tarfile=False):
    if not tarfile:
        print('install = {}'.format(compile_c_extension))
        if not libcpptraj_files:
            libcpptraj_files = try_updating_libcpptraj(
                    cpptraj_home=cpptraj_info.home_env,
                    compile_c_extension=compile_c_extension,
                    cpptraj_included=cpptraj_included,
                    openmp_flag=openmp_flag,
                    use_amberlib=use_amberlib
            )
    
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
            openmp_flag = ''
    
        extra_compile_args, extra_link_args = add_openmp_flag(openmp_flag,
            libcpptraj_has_openmp, extra_compile_args, extra_link_args)
    
        if sys.platform.startswith('linux'):
            # set rpath
            sys.stdout.write('set rpath to {}\n'.format(cpptraj_info.lib_dir))
            extra_link_args.append('-Wl,-rpath={}'.format(cpptraj_info.lib_dir))
            extra_compile_args.append('-Wl,-rpath={}'.format(cpptraj_info.lib_dir))
    
        check_cpptraj_version(cpptraj_info.include_dir, (4, 3, 1))
    
        pyxfiles, pxdfiles = get_pyx_pxd()

        if not is_released:
            from Cython.Build import cythonize
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
    
        library_dirs = [cpptraj_info.lib_dir, ]
    
        if sys.platform.startswith('darwin') and use_pip:
            # ship with libcpptraj.dylib in pytraj/lib/
            try:
                os.mkdir('pytraj/lib')
            except OSError:
                pass
    
            shutil.copy('{}/libcpptraj.dylib'.format(cpptraj_info.lib_dir), 'pytraj/lib')
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
                               include_dirs=[cpptraj_info.include_dir, pytraj_src],
                               extra_compile_args=extra_compile_args,
                               extra_link_args=extra_link_args)
            ext_modules.append(extmod)

        return ext_modules

def get_package_data(use_pip=False):
    sample_datafiles  = ["datafiles/ala3/Ala3.*",
                         "datafiles/tz2/tz2.*",
                         "datafiles/rna.pdb",
                         "datafiles/trpcage/trpcage*",
                         "datafiles/remd_ala2/*",
                         "datafiles/dpdp/DPDP*"]
    
    js_css_files = ['utils/progress-circle/css/*css',
                    'utils/progress-circle/*js',
                    'utils/css/*css',]
    
    _, pxdfiles = get_pyx_pxd()
    pkdata = pxdfiles + sample_datafiles + js_css_files
    
    if sys.platform.startswith('darwin') and use_pip:
        pkdata.append('lib/libcpptraj.dylib')
    return pkdata

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

if __name__ == '__main__':
    print(get_version_info())
