from __future__ import print_function, absolute_import
import os
import sys
from distutils.core import setup
from distutils.extension import Extension
from random import shuffle
from scripts.six import PY3 


# convert raw_input to input
PY3 = sys.version_info[0] == 3
if PY3:
    raw_input = input
# import/install Cython
try:
    import Cython.Distutils.build_ext
    from Cython.Build import cythonize
except:
    print ("There is no Cython")
    print ("try: pip install --upgrade git+git://github.com/cython/cython@master")
    try_cython = raw_input("install Cython? y/n ")

    if try_cython.upper() in ['Y', 'YES']:
        os.system("pip install --upgrade git+git://github.com/cython/cython@master")
    else:
        sys.exit("I can't install pytraj without cython")

# this setup.py file was adapted from setup.py file in Cython package
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

class PathError(Exception):
    def __init__(self, msg):
        pass

rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"

# find/install libcpptraj
try:
        cpptraj_dir = os.environ['CPPTRAJHOME'] 
        cpptraj_include = cpptraj_dir + "/src/"
        libdir = cpptraj_dir + "/lib/"
except:
    print ("have not set CPPTRAJHOME yet. \n")
    use_preshipped_lib = raw_input("use preshipped lib in ./Ambertools/V15.22b/cpptraj y/n \n")
    if use_preshipped_lib.upper() in ['Y', 'YES']:
        cpptraj_dir = rootname + "/Ambertools/V15.22b/cpptraj/"
        cpptraj_include = cpptraj_dir + "/src/"
        libdir = cpptraj_dir + "/lib"
    else:
        cpptraj_dir = ""

if not os.path.exists(cpptraj_dir):
    print ("cpptraj_dir does not exist")
    cpptraj_dir = raw_input("Please specify your cpptraj_dir: \n")
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = cpptraj_dir + "/lib/"
    #raise PathError("cpptraj_dir does not exist")

# get *.pyx files
pyxfiles = []
with open("PYXLIST.txt", 'r') as f:
    for line in f.readlines():
        if "#" not in line:
            pyxfiles.append(line.split("\n")[0])

#" use shuffle so we can use "python ./setup.py build" several times
# to make the compiling faster (really?)
shuffle(pyxfiles)

USE_PYX = True
if not USE_PYX:
    ext = ".cpp"
else:
    ext = ".pyx"

ext_modules = []
for ext_name in pyxfiles:
    pyxfile = pytraj_home + ext_name + ext

    # replace "/" by "." get module
    if "/" in ext_name:
        ext_name = ext_name.replace("/", ".")

    extmod = Extension("pytraj." + ext_name,
                    sources=[pyxfile],
                    libraries=['cpptraj'],
                    language='c++',
                    library_dirs=[libdir,],
                    include_dirs=[cpptraj_include, pytraj_home])
    ext_modules.append(extmod)

pxd_include_dirs = [
        directory for directory, dirs, files in os.walk('pytraj')
        if '__init__.pyx' in files or '__init__.pxd' in files
        or '__init__.py' in files
        ]

print (pxd_include_dirs)

pxd_include_patterns = [
        p+'/*.pxd' for p in pxd_include_dirs ] + [
        p+'/*.pyx' for p in pxd_include_dirs ]
         

setup_args = {}
# TODO : automatically discover packages
packages = [
        'pytraj',
        'pytraj.utils',
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
        ]

datalist = [[p[10:] for p in pxd_include_patterns]]

if __name__ == "__main__":
    setup(
        name="pytraj",
        version="0.1.beta",
        author="Hai Nguyen",
        author_email="hainm.comp@gmail.com",
        url="https://github.com/pytraj/pytraj",
        packages=packages,
        description="""Python API for cpptraj: a data analysis package for biomolecular simulation""",
        long_description=read("README.rst"),
        license = "BSD License",
        classifiers=[
                    'Development Status :: 5 - Production/Stable',
                    'Operating System :: Unix',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python 2.7',
                    'Programming Language :: Python 3.3',
                    'Programming Language :: Python 3.4',
                    'Programming Language :: Cython',
                    'Programming Language :: C',
                    'Programming Language :: C++',
                    'Topic :: Scientific/Engineering'],
        ext_modules = ext_modules,
        cmdclass = {'build_ext': Cython.Distutils.build_ext},
        package_data = {'pytraj' : ['data_sample/Ala3/Ala3.*',]},
    )

    print ()
    print ()
    print ("make sure to add libcpptraj to LD_LIBRARY_PATH before using pytraj")
    print ()
    print ()
