import os
from distutils.core import setup
from distutils.extension import Extension
import Cython.Distutils.build_ext
from Cython.Build import cythonize

# this setup.py file was adapted from setup.py file in Cython package
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"
cpptraj_dir = rootname + "/cpptraj/"
cpptraj_include = cpptraj_dir + "/src/"
libdir = cpptraj_dir + "/lib/"

#try:
#    try:
#        cpptraj_dir = cpptraj_include = os.environ['AMBERHOME'] + "/AmberTools/src/cpptraj/src/"
#        libdir = os.environ['AMBERHOME'] + "/lib"
#    except:
#        cpptraj_dir = cpptraj_include = os.environ['CPPTRAJHOME'] + "/src"
#        libdir = os.environ['CPPTRAJHOME'] + "/lib"
#except:
#    raise EnvironmentError("Must set AMBERHOME or CPPTRAJHOME")

# get *.pyx files
pyxfiles = []
with open("PYXLIST.txt", 'r') as f:
    for line in f.readlines():
        if "#" not in line:
            pyxfiles.append(line.split("\n")[0])

USE_PYX = True
if not USE_PYX:
    ext = ".cpp"
else:
    ext = ".pyx"

ext_modules = []
for ext_name in pyxfiles:
    pyxfile = pytraj_home + ext_name + ext

    # replace "/" by "." got get module
    if "/" in ext_name:
        ext_name = ext_name.replace("/", ".")

    extmod = Extension("pytraj." + ext_name,
                    [pyxfile],
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

pxd_include_patterns = [
        p+'/*.pxd' for p in pxd_include_dirs ] + [
        p+'/*.pyx' for p in pxd_include_dirs ]
         

setup_args = {}
packages = [
        'pytraj',
        'pytraj.utils',
        'pytraj.actions',
        'pytraj.analyses',
        ]

setup_args['package_data'] = {
        'pytraj' : [p[10:] for p in pxd_include_patterns],
        }

setup(
    name="PyCpptraj",
    version="0.1",
    author="Hai Nguyen",
    author_email="hainm.comp@gmail.com",
    url="https://github.com/hainm/pytraj",
    packages=packages,
    description="""Python wrapper for cpptraj""",
    long_description=read("../README.md"),
    license = "BSD License",
    platform = "",
    classifiers=[
                'Development Status :: Beta',
                'Operating System :: Linux',
                'Intended Audience :: Science/Research',
                'License :: BSD License',
                'Programming Language :: Python',
                'Programming Language :: Cython',
                'Programming Language :: C',
                'Programming Language :: C++',
                'Topic :: Scientific/Engineering'],
    ext_modules = ext_modules,
    cmdclass = {'build_ext': Cython.Distutils.build_ext},
    **setup_args
)
