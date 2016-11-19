#!/usr/bin/env python

"""build wheel for Linux, OSX
Requirement: python >= 3.4, conda, auditwheel and CPPTRAJHOME
Example: python build_wheel.py pytraj-1.0.9.tar.gz
Note: Recommend to build with centos 5 or 6 to maximize compatibility
"""
import os
import sys
import subprocess
from glob import glob

class PipBuilder(object):
    '''

    Parameters
    ----------
    tarfile : str, source file made by ./devtools/mkrelease
    pytraj_home : str, pytraj root folder
    python_versions : List[str], list of supported versions
        e.g.: ['2.7', '3.4', '3.5']

    Notes
    -----
    Update workflow in `run` method

    TODO
    ----
    Build wheel package from conda package?
    '''
    REQUIRED_PACKAGES = ['auditwheel']
    SUPPORTED_VERSIONS = ['2.7', '3.4', '3.5']
    MANY_LINUX_PYTHONS = dict(
       (py_version, '/opt/python/cp{py}-cp{py}m/bin/python'.format(py=py_version.replace('.', '')))
       for py_version in SUPPORTED_VERSIONS
    )

    def __init__(self,
                 tarfile,
                 pytraj_home,
                 python_versions,
                 use_manylinux=False,
                 cpptraj_dir=''):
        is_osx = sys.platform.startswith('darwin')
        self.python_versions = python_versions
        self.use_manylinux = use_manylinux
        if cpptraj_dir:
            self.cpptraj_dir = os.path.abspath(cpptraj_dir)
        else:
            self.cpptraj_dir = ''
        self.miniconda_root = self.find_miniconda_root()
        if self.use_manylinux:
            self.python_exe_paths = dict((py_version,
                PipBuilder.MANY_LINUX_PYTHONS[py_version]) 
                for py_version in self.python_versions)
        else:
            self.python_exe_paths = dict((py_version,
                                          '/'.join((self.miniconda_root,
                                                    '/envs/pytraj' + py_version,
                                                    'bin/python'))
                                          ) for py_version in self.python_versions)
        self.repair_exe = (pytraj_home + '/scripts/misc/fix_rpath_pip_osx.py' if is_osx else
                           'auditwheel repair')

    def run(self):
        self.check_cpptraj_and_required_libs()
        for python_version in self.python_versions:
            self.initialize_env(python_version)
            self.build_original_wheel(python_version)
            self.repair_wheel(python_version)
            self.validate_install(python_version)

    def initialize_env(self, python_version):
        if not self.use_manylinux:
            print('Using conda to create env')
            # use conda to create a new env
            env = 'pytraj' + python_version
            env_path = self.miniconda_root + '/envs/' + env
            if not os.path.exists(env_path):
                self.create_env(python_version)
            print('Building pytraj for Python {}'.format(python_version))
        else:
            print('Using python version from manylinux')

    def build_original_wheel(self, python_version):
        python = self.python_exe_paths[python_version]
        cmlist = '{python} -m pip wheel {tarfile}'.format(python=python, tarfile=tarfile).split()
        subprocess.check_call(cmlist)

    def find_miniconda_root(self):
        if self.use_manylinux:
            return ''
        else:
            command = "conda info | grep 'root environment'"
            output = subprocess.check_output(command, shell=True).decode()
            # e.g: outproot = "environment : /home/haichit/anaconda3  (writable)"
            return output.split()[3] + '/'

    def create_env(self, python_version):
        env = 'pytraj' + python_version
        sys.stdout.write('creating {} env'.format(env))
        cmlist = 'conda create -n {} python={} numpy nomkl --yes'.format(env, python_version)
        print(cmlist)
        subprocess.check_call(cmlist.split())
    
    def _get_wheel_file(self, python_version):
        return glob('*' + python_version.replace('.', '') +  '*.whl')[0]
    
    def repair_wheel(self, python_version):
        whl_name = self._get_wheel_file(python_version)
        command = '{} {}'.format(self.repair_exe, whl_name).split()
        subprocess.check_call(command)

    def validate_install(self, py_version):
        python_exe = self.python_exe_paths[py_version]
        env = 'pytraj' + py_version
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
        if self.use_manylinux:
            output = subprocess.check_output('{} -c "import pytraj as pt; print(pt)"'
                                           .format(python_exe, whl_file),
                                           shell=True)
            output = output.decode()
            print('Testing pytraj python={}'.format(py_version))
            print(output)
        else:
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
            print('PASSED: build test for {}'.format(whl_file))

    def check_cpptraj_and_required_libs(self):
        pytraj_home = os.path.abspath(os.path.dirname(__file__).strip('scripts'))
        cpptraj_home = os.environ.get("CPPTRAJHOME", '')
        if not cpptraj_home:
            print('CPPTRAJHOME is not set')
            ext = 'so' if not sys.platform.startswith('darwin') else 'dylib'
            if not self.cpptraj_dir:
                self.cpptraj_dir = pytraj_home + '/cpptraj/'
            suggested_libcpptraj = self.cpptraj_dir + '/lib/libcpptraj.' + ext
            print('Looking for {}'.format(suggested_libcpptraj))
            if os.path.exists(suggested_libcpptraj):
                print('Found')
                os.environ['CPPTRAJHOME'] = self.cpptraj_dir
                print('CPPTRAJHOME is set to {}'.format(self.cpptraj_dir))
            else:
                raise EnvironmentError('Must set CPPTRAJHOME')

        for package in self.REQUIRED_PACKAGES:
            try:
               __import__(package) 
            except ImportError:
                raise ImportError('require {}'.format(package))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Build wheel file to upload to pypi')
    parser.add_argument('tarfile')
    parser.add_argument('--py', default=None, help='Python version. Default: build all versions (2.7, 3.4, 3.5)')
    parser.add_argument('--cpptraj-dir', default='', help='cpptraj dir, optional')
    parser.add_argument('--manylinux-docker', action='store_true', help='If specified, use Python versions from manylinux')
    args = parser.parse_args()

    tarfile = os.path.abspath(args.tarfile)
    python_versions = PipBuilder.SUPPORTED_VERSIONS if args.py is None else [args.py,]
    # pytraj tar file
    pytraj_home = os.path.dirname(__file__).strip('scripts')
    builder = PipBuilder(tarfile=tarfile,
                         pytraj_home=pytraj_home,
                         python_versions=python_versions,
                         use_manylinux=args.manylinux_docker,
                         cpptraj_dir=args.cpptraj_dir)
    builder.run()
