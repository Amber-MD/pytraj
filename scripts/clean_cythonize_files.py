import os
from glob import glob

# get *.pyx files
pyxfiles = []
pxd_include_dirs = [
    directory for directory, dirs, files in os.walk('pytraj')
    if '__init__.pyx' in files or '__init__.pxd' in files
    or '__init__.py' in files
]

pxd_include_patterns = [
    p + '/*.pxd' for p in pxd_include_dirs]

for p in pxd_include_dirs:
    pyxfiles.extend([ext.split(".")[0] for ext in glob(p + '/*.pyx') if '.pyx' in ext])

print("move old cythonized files to ./trash folder")
for ext_name in pyxfiles:
    pyxfile = ext_name + ".cpp"
    sofile = ext_name + "*.so"

    for f in [pyxfile, sofile]:
        do_this = "mv %s trash" % f
        os.system(do_this)
