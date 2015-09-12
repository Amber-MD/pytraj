import sys
import subprocess
from subprocess import CalledProcessError
from glob import glob
from itertools import chain
# import ./scripts
from scripts import setup_for_amber

if sys.version_info < (2, 6):
    sys.stderr.write('You must have at least Python 2.6 for pytraj\n')
    sys.exit(0)

import os
from distutils.core import setup, Command
from distutils import ccompiler
from distutils.extension import Extension
from random import shuffle
import time
from time import sleep

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

pytraj_version = read("pytraj/__version__.py").split("=")[-1]
pytraj_version = pytraj_version.replace('"', '', 10)
rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"

if len(sys.argv) == 2 and sys.argv[1] == 'install':
    do_install = True
else:
    do_install = False

if len(sys.argv) == 2 and sys.argv[1] == 'build':
    do_build = True
else:
    do_build= False


# check/install Cython
cmdclass = {}
try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
    has_cython = True
    cmdclass['build_ext'] = build_ext
except ImportError:
    has_cython = False
    from distutils.command.build_ext import build_ext 
    cmdclass['build_ext'] = build_ext
    #sys.stderr.write('You must have Cython installed to install pytraj\n')
    #sys.exit(0)

# check AMBERHOME
try:
    amberhome = os.environ['AMBERHOME']
    has_amberhome = True
except KeyError:
    has_amberhome = False

# check CPPTRAJHOME or "./cpptraj" folder
try:
    cpptrajhome = os.environ['CPPTRAJHOME']
    has_cpptrajhome = True
except KeyError:
    has_cpptrajhome = False

has_cpptraj_in_current_folder = os.path.exists("./cpptraj/")

if has_cpptrajhome:
    # use libcpptraj and header files in CPPTRAJHOME (/lib, /src)
    cpptraj_dir = cpptrajhome
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = cpptrajhome + "/lib/"
elif has_cpptraj_in_current_folder:
    cpptraj_dir = os.path.abspath("./cpptraj/")
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = cpptraj_dir + "/lib/"
# turn off using AMBERHOME since pytraj uses cpptraj-dev
#elif has_amberhome:
#    # use libcpptraj and header files in AMBERHOME
#    cpptraj_dir = amberhome + "/AmberTools/src/cpptraj/"
#    cpptraj_include = cpptraj_dir + "/src/"
#    libdir = amberhome + "/lib/"
else:
    # use libcpptraj and header files in PYTRAJHOME
    # ./cpptraj/src
    # ./cpptraj/lib

    nice_message = """
    Trying to dowload and build libcpptraj for you. (5-10 minutes)
    (check ./cpptraj/ folder after installation)

    To avoid auto-installation
    --------------------------
    Must set CPPTRAJHOME or installing ./cpptraj/ in current folder.

    If you want to manually install `libcpptraj`, you can download cpptraj
    development version from here: https://github.com/Amber-MD/cpptraj

    (    git clone https://github.com/Amber-MD/cpptraj/
         cd cpptraj
         export CPPTRAJHOME=`pwd`
         ./configure -shared gnu
         make libcpptraj    )

    and then go back to pytraj folder:
    python ./setup.py install

    ...
    """
    if do_install or do_build:
        print (nice_message)
        sleep(3)
        try:
            subprocess.check_call(['sh', './installs/install_cpptraj_git.sh'])
        except CalledProcessError:
            sys.stderr.write('can not install libcpptraj, you need to install it manually \n')
            sys.exit(1)
    cpptraj_dir = os.path.join(rootname, "cpptraj")
    cpptraj_include = os.path.join(cpptraj_dir, 'src')
    libdir =  os.path.join(cpptraj_dir, 'lib')

# get *.pyx files
pxd_include_dirs = [
        directory for directory, dirs, files in os.walk('pytraj')
        if '__init__.pyx' in files or '__init__.pxd' in files
        or '__init__.py' in files
        ]

pxd_include_patterns = [ 
        p + '/*.pxd' for p in pxd_include_dirs ]

pyxfiles = []
for p in pxd_include_dirs:
    pyxfiles.extend([ext.split(".")[0] for ext in glob(p + '/*.pyx') if '.pyx' in ext])

# check command line
extra_compile_args=['-O0',]
extra_link_args=['-O0',]

list_of_libcpptraj = glob(os.path.join(libdir, 'libcpptraj') + '*')
if not list_of_libcpptraj:
    if do_install or do_build:
        raise RuntimeError('can not find libcpptraj in $CPPTRAJHOME/lib or ./cpptraj/lib folder. '
                           'You need to install ``libcpptraj`` manually. '
                           )

openmp_str = "openmp"
if openmp_str in sys.argv:
    # python ./setup.py build openmp
    # make sure to update Makefile in $AMBERHOME/AmberTools/src
    # if changing '-openmp' to something else
    with_openmp = True
    sys.argv.remove(openmp_str)
else:
    with_openmp = False 

if with_openmp:
    extra_compile_args.append("-fopenmp")
    extra_link_args.append("-fopenmp")

KeyErrorTXT = """
Can not use -faster_build with `install`,
try  "python ./setup.py build faster_build
then "python ./setup.py install" 
"""

faster_build_str = "faster"
if faster_build_str in sys.argv:
    # try using multiple cores
    faster_build = True
    sys.argv.remove(faster_build_str)
    if "install" in sys.argv:
        sys.stderr.write(KeyErrorTXT)
        sys.exit(0)
else:
    faster_build = False

# since we added "INSTALLTYPE" after setup.py file, we need
# to remove it if having one
installtype = os.environ.get("INSTALLTYPE", "")
try:
    sys.argv.remove(installtype)
except ValueError:
    pass

ext_modules = []
for ext_name in pyxfiles:
    if has_cython:
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
                    library_dirs=[libdir,],
                    include_dirs=[cpptraj_include, pytraj_home],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args)

    extmod.cython_directives = {
            'embedsignature':True,
            'boundscheck': False,
            'wraparound':False,
            }
    ext_modules.append(extmod)

#shuffle(ext_modules)

setup_args = {}
packages = [
        'pytraj',
        'pytraj.utils',
        'pytraj.actions',
        'pytraj.analyses',
        'pytraj.datasets',
        'pytraj.externals',
        'pytraj.trajs',
        'pytraj.datafiles',
        'pytraj.datafiles.Ala3',
        'pytraj.datafiles.tz2',
        'pytraj.plot',
        'pytraj.math',
        'pytraj.core',
        'pytraj.parallel',
        'pytraj.cluster',
        'pytraj.sandbox',
        ]

pylen = len('pytraj') + 1
pxdlist = [p.replace("pytraj/", "") for p in pxd_include_patterns]
sample_data = ["datafiles/Ala3/Ala3.*", "datafiles/tz2/tz2.*"]
datalist = pxdlist +  sample_data

# compare to "setup_for_amber" script
package_match = (sorted(packages) == sorted(setup_for_amber.packages))
datalist_match = (sorted(datalist) == sorted(setup_for_amber.datalist))

if not package_match:
    sys.stderr.write("packages mistmatch. Make sure to update ./scripts/setup_for_amber.py\n")
    sys.exit(1)

if not datalist_match:
    sys.stderr.write("datalist mistmatch\n")
    sys.exit(1)
    
def build_func(my_ext):
    return setup(name="pytraj",
        version=pytraj_version,
        author="Hai Nguyen",
        author_email="hainm.comp@gmail.com",
        url="https://github.com/Amber-MD/pytraj",
        packages=packages,
        description="""Python API for cpptraj: a data analysis package for biomolecular simulation""",
        license = "BSD License",
        classifiers=[
                    'Development Status :: 4 - Beta',
                    'Operating System :: Unix',
                    'Operating System :: MacOS',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 3.3',
                    'Programming Language :: Python :: 3.4',
                    'Programming Language :: Cython',
                    'Programming Language :: C',
                    'Programming Language :: C++',
                    'Topic :: Scientific/Engineering'],
        ext_modules = my_ext,
        package_data = {'pytraj' : datalist},
        cmdclass = cmdclass,
        )

def remind_ld_lib_path(build_tag, libdir):
    if build_tag:
        from scripts.acsii_art import batman
        print ("")
        print ("")
        print (batman)
        libdir = os.path.abspath(libdir)
        print ("make sure to add %s to your LD_LIBRARY_PATH" % libdir)
        print ("")
    else:
        print ("not able to install pytraj")

if __name__ == "__main__":
    if not faster_build:
        build_tag = build_func(ext_modules)
        if do_install:
            remind_ld_lib_path(build_tag, libdir)
    else:
        from multiprocessing import cpu_count
        n_cpus = cpu_count()
        print('number of available cpus = %s' % n_cpus)
        num_each = int(len(ext_modules)/n_cpus)

        sub_ext_modules_list = []
        # there is idiom to do this but I am too lazy to think
        for i in range(n_cpus):
            if i != n_cpus-1:
                sub_ext_modules_list.append(ext_modules[i*num_each:num_each*(i+1)])
            else:
                sub_ext_modules_list.append(ext_modules[num_each*i:])

        from multiprocessing import Pool
        pool = Pool(n_cpus)
        pool.map(build_func, sub_ext_modules_list)
