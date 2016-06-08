import os
import sys
import Cython
from Cython.Build import cythonize
from glob import glob

if Cython.__version__ < '0.21':
    sys.stderr.write('require cython version >=0.21')
    sys.exit(1)

DEBUG = False

pxd_include_dirs = [
    directory for directory, dirs, files in os.walk('pytraj') if '__init__.pyx'
    in files or '__init__.pxd' in files or '__init__.py' in files
]

pyxfiles = []
for p in pxd_include_dirs:
    pyxfiles.extend([ext.split(".")[0] for ext in glob(p + '/*.pyx')
                     if '.pyx' in ext])

cython_directives = {
    'embedsignature': True,
    'boundscheck': False,
    'wraparound': False,
}

if DEBUG:
    cython_directives.update({
        'profile': True,
        'linetrace': True,
        'binding': True})
    define_macros = [('CYTHON_TRACE', 1), ]
else:
    define_macros = []

cythonize(
    [pfile + '.pyx' for pfile in pyxfiles],
    nthreads=int(os.environ.get('NUM_THREADS', 4)),
    compiler_directives=cython_directives,
)
