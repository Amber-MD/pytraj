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
    #sys.exit(0)

rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"
amber_home = os.environ['AMBERHOME']

cpptraj_dir = amber_home + "/AmberTools/src/cpptraj/"
cpptraj_include = cpptraj_dir + "/src/"
libdir = amber_home + "/lib/"

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

# make random list so we can run many `python setup.py` at the same times
# TODO: need to compile parallely.
shuffle(pyxfiles)

extra_compile_args=['-O0', '-ggdb']
extra_link_args=['-O0', '-ggdb']

if "-openmp" in sys.argv:
    with_openmp = True
    sys.argv.remove("-openmp")
else:
    with_openmp = False 

if with_openmp:
    extra_compile_args.append("-fopenmp")
    extra_link_args.append("-fopenmp")

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
                    library_dirs=[libdir,],
                    include_dirs=[cpptraj_include, pytraj_home],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args)

    extmod.cython_directives = {
            'embedsignature':True,
            'boundscheck': False,
            }
    ext_modules.append(extmod)

pxd_include_dirs = [
        directory for directory, dirs, files in os.walk('pytraj')
        if '__init__.pyx' in files or '__init__.pxd' in files
        or '__init__.py' in files
        ]

pxd_include_patterns = [
        p+'/*.pxd' for p in pxd_include_dirs ]

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
sample_data = ["data_sample/Ala3/Ala3.*",]
html_data = ["html/static/*"] 
datalist = pxdlist +  sample_data + html_data

if __name__ == "__main__":
    setup(name="pytraj",
        version="0.1.2.dev0",
        author="Hai Nguyen",
        author_email="hainm.comp@gmail.com",
        url="https://github.com/pytraj/pytraj",
        packages=packages,
        description="""Python API for cpptraj: a data analysis package for biomolecular simulation""",
        long_description=read("../README.rst"),
        license = "BSD License",
        classifiers=[
                    'Development Status :: 5 - Production/Stable',
                    'Operating System :: Unix',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 2',
                    'Programming Language :: Python :: 3',
                    'Programming Language :: Cython',
                    'Programming Language :: C',
                    'Programming Language :: C++',
                    'Topic :: Scientific/Engineering'],
        ext_modules = ext_modules,
        package_data = {'pytraj' : datalist},
        cmdclass = cmdclass,
    )
