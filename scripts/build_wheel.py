#!/usr/bin/env python

"""build wheel for Linux, OSX
Requirement: python >= 3.4, conda, auditwheel and CPPTRAJHOME
Example: python build_wheel.py pytraj-1.0.9.tar.gz
Note: Recommend to build with centos 5 or 6 to maximize compatibility
"""
import os, sys
import subprocess
from glob import glob

REQUIRED_PACKAGES = ['auditwheel']
SUPPORTED_VERSION = ['2.7', '3.4', '3.5']
assert os.environ.get("CPPTRAJHOME"), "must set CPPTRAJHOME"

def check_package():
    for package in REQUIRED_PACKAGES:
        try:
           __import__(package) 
        except ImportError:
            raise ImportError('require {}'.format(package))

def find_miniconda_root():
    command = "conda info | grep 'root environment'"
    output = subprocess.check_output(command, shell=True).decode()
    # e.g: outproot = "environment : /home/haichit/anaconda3  (writable)"
    return output.split()[3] + '/'

def build(tarfile,
          pytraj_home,
          miniconda_root,
          osx_rpath_script,
          python_versions):
    working_version = '.'.join((str(i) for i in sys.version_info[:2]))
    print("working_version = {}".format(working_version))
    for python_version in python_versions:
        print('building pytraj for Python {}'.format(python_version))
        env = 'pytraj' + python_version

        # create whl file for each python version
        print(miniconda_root + env)
        if not os.path.exists(miniconda_root + env):
            sys.stdout.write('creating {} env'.format(env))
            cmlist = 'conda create -n {} python={} numpy nomkl --yes'.format(env, python_version)
            print(cmlist)
            subprocess.check_call(cmlist.split())
    
        python = '/'.join((miniconda_root, env, 'bin/python'))
    
        if not (sys.platform.startswith('darwin') and sys.version_info == (2, 7)):
            cmlist = '{python} -m pip wheel {tarfile}'.format(python=python, tarfile=tarfile).split()
            subprocess.check_call(cmlist)

        whl_name = glob('*' + python_version.replace('.', '') +  '*.whl')[0]
    
        if sys.platform.startswith('darwin') and sys.version_info > (2, 7):
            # only have fix for python >= 3.4
            exe = osx_rpath_script
        elif sys.platform.startswith('linux'):
            exe = 'auditwheel repair'
        else:
            exe = 'echo: ignore '
    
        cmlist = '{} {}'.format(exe, whl_name).split()
        subprocess.check_call(cmlist)
        test_install(python, python_version, env)

def test_install(python_exe, py_version, env):
    # e.g: change 2.7 to 27
    print('Testing pytraj build')
    version = py_version.replace('.', '')
    cwd = os.getcwd()
    pattern = cwd + '/wheelhouse/pytraj-*-cp{}-*.whl'.format(version)
    whl_file = glob(pattern)[0]
    try:
        subprocess.check_call('{} -m pip uninstall pytraj -y'.format(python_exe).split())
    except subprocess.CalledProcessError:
        pass
    subprocess.check_call('{} -m pip install {}'.format(python_exe, whl_file).split())
    try:
        print('Installing numpy')
        subprocess.check_call('{} -c "import numpy"'.format(python_exe), shell=True)
    except subprocess.CalledProcessError:
        subprocess.check_call('conda install numpy nomkl -y -n {}'.format(env),
                              shell=True)
    output = subprocess.check_output('{} -c "import pytraj as pt; print(pt)"'
                                   .format(python_exe, whl_file),
                                   shell=True)
    output = output.decode()
    expected_output = 'envs/{env}/lib/python{py_version}/site-packages/pytraj'.format(env=env,
            py_version=py_version)
    assert expected_output in output
    print('PASSED: build test')

if __name__ == '__main__':
    check_package()
    import argparse
    parser = argparse.ArgumentParser('Build wheel file to upload to pypi')
    parser.add_argument('tarfile')
    parser.add_argument('--py', default=None, help='Python version. Default: build all versions (2.7, 3.4, 3.5)')
    args = parser.parse_args()
    tarfile = os.path.abspath(args.tarfile)
    python_versions = SUPPORTED_VERSION if args.py is None else [args.py,]

    # pytraj tar file
    miniconda_root = find_miniconda_root() + '/envs/'
    pytraj_home = os.path.dirname(__file__).strip('scripts')
    osx_rpath_script = pytraj_home + '/scripts/misc/fix_rpath_pip_osx.py'
    build(tarfile=tarfile,
          pytraj_home=pytraj_home,
          miniconda_root=miniconda_root,
          osx_rpath_script=osx_rpath_script,
          python_versions=python_versions)
