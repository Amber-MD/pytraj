# TODO: make clean setup file
from __future__ import print_function, absolute_import
import os
import sys
from distutils.core import setup
from distutils import ccompiler
from distutils.extension import Extension
from random import shuffle
import time

# global
cpptraj_include = "" 
cpptraj_dir = "" 
libdir = ""

def find_amberhome():
    try:
        amber_home = os.environ['AMBERHOME']
        return amber_home
    except:
        return None

def find_libnetcdef():
    compiler=ccompiler.new_compiler()
    _lib_dirs = os.environ['PATH'].split(":") 
    home_dir = os.environ['HOME']
    lib_dirs = _lib_dirs + [src + "/lib" for src in _lib_dirs ] +  [home_dir + "/anaconda/lib/", ]
    libnetcdf = compiler.find_library_file(lib_dirs, 'netcdf')
    return libnetcdf

has_netcdf = True if find_libnetcdef() else False

def install_cpptraj_git():
    global cpptraj_include
    global cpptraj_dir
    global libdir

    print ("download from http://github.com/mojyt/cpptraj")
    cpptraj_dir = rootname + "/cpptraj/"
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = cpptraj_dir + "/lib"
    old_file = "install_cpptraj_git.sh"
    new_file = "_" + old_file
    add_netcdf_to_install_file(old_file, new_file, has_netcdf)
    os.system("sh ./installs/" + new_file)

def add_netcdf_to_install_file(fh1, fh2, has_netcdf):
    # assume libnetcdf*so in $NETCDF_HOME/lib
    # and header file is in $NETCDF_HOME/include
    fh1 = "./installs/" + fh1
    fh2 = "./installs/" + fh2
    # always use $AMBERHOME/lib for libnetcdf
    # if not, use others
    if find_amberhome() is None:
        if find_libnetcdef():
            libnetcdf_dir = find_libnetcdef().split("lib")[0]
        else:
            libnetcdf_dir = None
    else:
        libnetcdf_dir = find_amberhome()

    # add netcdf flag
    if has_netcdf:
        new_chunk = "-shared --with-netcdf=%s gnu" % libnetcdf_dir
    else:
        new_chunk = "-shared -nonetcdf gnu"
    with open(fh1, 'r') as _fh1, open(fh2, 'w') as _fh2:
        txt = _fh1.read()
        txt = txt.replace("-shared gnu", new_chunk)
        _fh2.write(txt)

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# convert raw_input to input
PY3 = sys.version_info[0] == 3
if PY3:
    raw_input = input

# check Cython
try:
    import Cython.Distutils.build_ext
    from Cython.Build import cythonize
    has_cython = True
except:
    has_cython = False

rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"

# find/install libcpptraj
try:
    cpptraj_dir = os.environ['CPPTRAJHOME'] 
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = cpptraj_dir + "/lib/"
    has_cpptrajhome = True
except:
    has_cpptrajhome = False
    using_pip = "pip" in os.path.dirname(__file__)
    print (os.path.basename(os.path.dirname(__file__)))
    print ("using_pip = %s" % using_pip)
    if using_pip:
        sys.exit()
        os.system("pip install -d . pytraj")
        os.chdir("./pytraj")
        os.system("python ./setup.py install")
    else:
        print ()
        print ("You have not yet set CPPTRAJHOME. \n")
        print ("To avoid below message everytime you install/build ..., just set CPPTRAJHOME")
        print ("you have two options: \n")
        print ()
        print ("Option 1. Quit and export/setenv CPPTRAJHOME and/or install libcpptraj")
        print ("        + export CPPTRAJHOME=your_favorite_dir (if using bash) / setenv CPPTRAJHOME your_favorite_dir")
        print ("        + if not install libcpptraj yet, please follow")
        print ("        + (cd $CPPTRAJHOME)")
        print ("        + (./configure -shared gnu)")
        print ("        + (make libcpptraj)\n")
        print ("Option 2. use cpptraj development version in github: https://github.com/mojyt/cpptraj")
        print ("         (Bonus: I will automatically install it for you, don't worry)\n")
        print ("choose Quit (1) or Stay (2)? \n")
        
        answer = raw_input("\n")
        if answer.lower() in ['q', 'quite', '1']:
            print ("Bye bye ^_^")
            print ("Quit ....")
            from scripts.acsii_art import batman
            print (batman)
            sys.exit()
        else:
            use_cpptraj_git = True

        if use_cpptraj_git:
            install_cpptraj_git()

if not os.path.exists(cpptraj_dir):
    print ("cpptraj_dir does not exist")
    cpptraj_dir = raw_input("Please specify your cpptraj_dir: \n")
    cpptraj_include = cpptraj_dir + "/src/"
    libdir = cpptraj_dir + "/lib/"

# get *.pyx files
pyxfiles = []
with open("pyxlist.txt", 'r') as f:
    for line in f.readlines():
        if "#" not in line:
            pyxfiles.append(line.split("\n")[0])

#" use shuffle so we can use "python ./setup.py build" several times
# to make the compiling faster (really?)
shuffle(pyxfiles)

ext = ".pyx"
ext_modules = []
for ext_name in pyxfiles:
    if has_cython:
        pyxfile = pytraj_home + ext_name + ext
    else:
        pyxfile = pytraj_home + ext_name + ".cpp"

    # replace "/" by "." get module
    if "/" in ext_name:
        ext_name = ext_name.replace("/", ".")

    sources = [pyxfile]
    #if "Action_Mask_2" in ext_name:
    #    sources.append("./pytraj/src/Action_Mask_2.cpp")

    extmod = Extension("pytraj." + ext_name,
                    sources=sources,
                    libraries=['cpptraj'],
                    language='c++',
                    library_dirs=[libdir,],
                    include_dirs=[cpptraj_include, pytraj_home],
                    extra_compile_args=['-O0', '-ggdb'],
                    extra_link_args=['-O0', '-ggdb'])

    extmod.cython_directives = {
            'embedsignature':True,
            'boundscheck': False,
            #'profile': True,
            #'linetrace' : True
            }
    ext_modules.append(extmod)

pxd_include_dirs = [
        directory for directory, dirs, files in os.walk('pytraj')
        if '__init__.pyx' in files or '__init__.pxd' in files
        or '__init__.py' in files
        ]

#print (pxd_include_dirs)

pxd_include_patterns = [
        p+'/*.pxd' for p in pxd_include_dirs ] + [
        p+'/*.pyx' for p in pxd_include_dirs ]
         

setup_args = {}
# TODO : automatically discover packages
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
        ]

pylen = len('pytraj') + 1
#print (pylen)
# remove `pytraj` name 
#datalist = [[p[pylen:] for p in pxd_include_patterns]]
datalist = [p for p in pxd_include_patterns]
#print (datalist)
sample_data = ["data_sample/Ala3/Ala3.*",]
html_data = ["html/static/*"] 
datalist = datalist +  sample_data + html_data
#print (datalist)

if __name__ == "__main__":
    setup(
        name="pytraj",
        version="0.1.0.2pre",
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
                    'Programming Language :: Python :: 2',
                    'Programming Language :: Python :: 3',
                    'Programming Language :: Cython',
                    'Programming Language :: C',
                    'Programming Language :: C++',
                    'Topic :: Scientific/Engineering'],
        ext_modules = ext_modules,
        cmdclass = {'build_ext': Cython.Distutils.build_ext},
        package_data = {'pytraj' : datalist},
    )

    print ()
    print ()
    print ("make sure to add libcpptraj to LD_LIBRARY_PATH before using pytraj")
    print (libdir)
    print ()
    print ("or move libcpptraj.so to existing LD_LIBRARY_PATH")
    print ()
    print ()
