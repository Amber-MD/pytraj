#!/usr/bin/env python

"""build wheel for Linux, OSX

Requirement: conda, auditwheel CPPTRAJHOME
"""
import os, sys
from glob import glob

assert os.environ.get("CPPTRAJHOME"), "must set CPPTRAJHOME"

# pytraj tar file
miniconda_root = sys.exec_prefix  + '/envs/'
pytraj_home = "/Users/haichit/pytraj_github/pytraj/"
fix_osx = pytraj_home + '/scripts/misc/fix_rpath_pip_osx.py'

pytraj_envs = ['pytraj2.7', 'pytraj3.4', 'pytraj3.5']
# pytraj_envs = ['pytraj3.4', 'pytraj3.5']

def build(pytar,
          miniconda_root=miniconda_root,
          pytraj_home=pytraj_home,
          fix_osx=fix_osx,
          pytraj_envs=pytraj_envs):
    for pyenv in pytraj_envs:
        pyversion = pyenv.strip('pytraj')
        if not os.path.exists(miniconda_root + pyenv):
            sys.stdout.write('creating {} env'.format(pyenv))
            os.system('conda create -n {} python={} --yes'.format(pyenv, pyversion))
    
        python = '/'.join((miniconda_root, pyenv, 'bin/python'))
    
        if not (sys.platform.startswith('darwin') and sys.version_info == (2, 7)):
            os.system('{python} -m pip wheel {pytar}'.format(python=python, pytar=pytar))
        whl_name = glob('*' + pyversion.replace('.', '') +  '*.whl')[0]
    
        if sys.platform.startswith('darwin') and sys.version_info > (2, 7):
            # only have fix for python >= 3.4
            exe = fix_osx
        elif sys.platform.startswith('linux'):
            exe = 'auditwheel repair'
        else:
            exe = 'echo: ignore '
    
        os.system('{} {}'.format(exe, whl_name))

if __name__ == '__main__':
    pytar = sys.argv[1]
    pytraj_home = os.path.abspath(__file__).replace('scripts/build_wheel.py', '')
    fix_osx = pytraj_home + '/scripts/misc/fix_rpath_pip_osx.py'
    build(pytar=pytar, pytraj_home=pytraj_home, fix_osx=fix_osx)

