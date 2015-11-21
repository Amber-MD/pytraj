import sys
from glob import glob

if sys.version_info < (2, 6):
    sys.stderr.write('You must have at least Python 2.6 for pytraj\n')
    sys.exit(0)

import os
from distutils.core import setup, Command
from distutils import ccompiler
from distutils.extension import Extension
from random import shuffle
import time


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

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
    # sys.exit(0)

rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"

try:
    amberhome = os.environ['AMBERHOME']
    has_amberhome = True
except:
    has_amberhome = False

try:
    cpptrajhome = os.environ['CPPTRAJHOME']
    has_cpptrajhome = True
except:
    has_cpptrajhome = False

if has_amberhome:
    cpptraj_dir = amberhome + "/AmberTools/src/cpptraj/"
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = amberhome + "/lib/"
elif has_cpptrajhome:
    cpptraj_dir = cpptrajhome
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = cpptrajhome + "/lib/"
else:
    cpptraj_dir, cpptraj_include, libdir = None, None, None
    nice_message = """
    Must set AMBERHOME or CPPTRAJHOME. If both AMBERHOME and CPPTRAJHOME are set,
    pytraj will give priority to AMBERHOME.

    If you don't have AmberTools or cpptraj, you can install pytraj by
    one of two ways here:

    1. Download AmberTools15 (or later version)
    First, get a free version from: http://ambermd.org/#AmberTools
    then:
        cd $AMBERHOME/AmberTools/src/
        make pytraj

    2. if you just want a standalone cpptraj version, you can download
    development version from here: https://github.com/mojyt/cpptraj

    (    git clone https://github.com/mojyt/cpptraj
         cd cpptraj
         export CPPTRAJHOME=`pwd`
         ./configure -shared gnu
         make libcpptraj    )

    and then go back to pytraj folder:
    python ./setup.py install
    """
    sys.stderr.write(nice_message)
    sys.exit(0)

# find and give warning if not having libcpptraj?

# get *.pyx files
pyxfiles = []
f = open('pyxlist.txt', 'r')
try:
    for line in f.readlines():
        if "#" not in line:
            pyxfiles.append(line.split("\n")[0])
finally:
    f.close()

ext_modules = []
for ext_name in pyxfiles:
    if has_cython:
        ext = ".pyx"
    else:
        ext = ".cpp"
    pyxfile = pytraj_home + ext_name + ext

    # replace "/" by "." get module
    if "/" in ext_name:
        ext_name = ext_name.replace("/", ".")

    sources = [pyxfile]

    extmod = Extension("pytraj." + ext_name,
                       sources=sources,
                       libraries=['cpptraj'],
                       language='c++',
                       library_dirs=[libdir, ],
                       include_dirs=[cpptraj_include, pytraj_home],
                       extra_compile_args=['-O0', '-ggdb', ],
                       extra_link_args=['-O0', '-ggdb', ],)

    extmod.cython_directives = {
        'embedsignature': True,
        'boundscheck': False,
    }
    ext_modules.append(extmod)

pxd_include_dirs = [
    directory for directory, dirs, files in os.walk('pytraj')
    if '__init__.pyx' in files or '__init__.pxd' in files
    or '__init__.py' in files
]

pxd_include_patterns = [
    p + '/*.pxd' for p in pxd_include_dirs]

setup_args = {}
packages = [
    'pytraj',
    'pytraj.utils',
    'pytraj.html',
    'pytraj.actions',
    'pytraj.analyses',
    'pytraj.datasets',
    'pytraj.externals',
    'pytraj.parms',
    'pytraj.clusters',
    'pytraj.trajs',
    'pytraj.gdt',
    'pytraj.data_sample',
    'pytraj.data_sample.Ala3',
    'pytraj.plots',
]

pylen = len('pytraj') + 1
pxdlist = [p.replace("pytraj/", "") for p in pxd_include_patterns]
sample_data = ["data_sample/Ala3/Ala3.*", ]
html_data = ["html/static/*"]
datalist = pxdlist + sample_data + html_data


def build_func(my_ext):
    setup(name="pytraj",
          version="0.1.2.dev1",
          author="Hai Nguyen",
          author_email="hainm.comp@gmail.com",
          url="https://github.com/pytraj/pytraj",
          packages=packages,
          description="""Python API for cpptraj: a data analysis package for biomolecular simulation""",
          long_description=read("../README.rst"),
          license="BSD License",
          classifiers=[
              'Development Status :: 4 - Beta',
              'Operating System :: Unix',
              'Operating System :: MacOS',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Programming Language :: Python :: 2.6'
              'Programming Language :: Python :: 2.7',
              'Programming Language :: Python :: 3.3',
              'Programming Language :: Python :: 3.4',
              'Programming Language :: Cython',
              'Programming Language :: C',
              'Programming Language :: C++',
              'Topic :: Scientific/Engineering'],
          ext_modules=my_ext,
          package_data={'pytraj': datalist},
          cmdclass=cmdclass,
          )


def test_func(my_ext):
    print(len(my_ext))

if __name__ == "__main__":
    from multiprocessing import cpu_count
    n_cpus = cpu_count()
    num_each = int(len(ext_modules) / n_cpus)

    sub_ext_modules_list = []
    for i in range(n_cpus):
        if i != n_cpus - 1:
            sub_ext_modules_list.append(ext_modules[i * num_each:num_each * (i + 1)])
        else:
            sub_ext_modules_list.append(ext_modules[num_each * i:])

    from multiprocessing import Pool
    pool = Pool(n_cpus)
    pool.map(build_func, sub_ext_modules_list)
