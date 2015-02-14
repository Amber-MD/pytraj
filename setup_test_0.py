from __future__ import print_function, absolute_import
import os
import sys
from glob import glob
from distutils.core import setup
from distutils.extension import Extension
from distutils import ccompiler
from random import shuffle

# this setup.py file was adapted from setup.py file in Cython package
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def find_netcdef():
    compiler=ccompiler.new_compiler()
    _lib_dirs = os.environ['PATH'].split(":") 
    home_dir = os.environ['HOME']
    lib_dirs = _lib_dirs + [src + "/lib" for src in _lib_dirs ] +  [home_dir + "/anaconda3/lib/", ]
    libnetcdf_dir = compiler.find_library_file(lib_dirs, 'netcdf')
    return libnetcdf_dir

print (find_netcdef())

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
    print ("use cpptraj git: git clone http://github.com/mojyt/cpptraj?")
    use_cpptraj_git = raw_input("y/n \n")
    if use_cpptraj_git.upper() in ['Y', 'YES']:
        print ("download from http://github.com/mojyt/cpptraj")
        cpptraj_dir = rootname + "/cpptraj/"
        cpptraj_include = cpptraj_dir + "/src/"
        libdir = cpptraj_dir + "/lib"
        os.system("./installs/install_cpptraj_git.sh")
    else:
        use_preshipped_lib = raw_input("use preshipped lib in ./Ambertools/dev/cpptraj y/n \n")
        if use_preshipped_lib.upper() in ['Y', 'YES']:
            cpptraj_dir = rootname + "/Ambertools/dev/cpptraj/"
            cpptraj_include = cpptraj_dir + "/src/"
            libdir = cpptraj_dir + "/lib"
            os.system("sh ./installs/install_cpptraj.sh")
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
        version="0.1.beta.7",
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
                    'Programming Language :: Python',
                    'Programming Language :: Python',
                    'Programming Language :: Python',
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
    print (libdir)
    print ()
    print ()
