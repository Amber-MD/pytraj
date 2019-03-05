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

sys.path.insert(0, os.path.abspath(__file__))
from utils import temporarily_move_libcpptraj


class PipBuilder(object):
    '''

    Parameters
    ----------
    tarfile : str, source file made by ./devtools/mkrelease
    pytraj_home : str, pytraj root folder
    python_versions : List[str], list of supported versions
        e.g.: ['2.7', '3.4', '3.5', '3.6']

    Notes
    -----
    Update workflow in `run` method

    TODO
    ----
    Build wheel package from conda package?
    '''
    REQUIRED_PACKAGES = ['auditwheel']
    SUPPORTED_VERSIONS = ['2.7mu', '2.7', '3.4', '3.5', '3.6']

    if sys.platform.startswith('darwin'):
        REQUIRED_PACKAGES.append('conda_build')
        SUPPORTED_VERSIONS.remove('2.7mu')

    MANY_LINUX_PYTHONS = dict((py_version,
                               '/opt/python/cp{py}-cp{py}m/bin/python'.format(
                                   py=py_version.replace('.', '')))
                              for py_version in SUPPORTED_VERSIONS)
    # wide-unicode (to be compatible with miniconda/anaconda python)
    if '2.7mu' in MANY_LINUX_PYTHONS:
        MANY_LINUX_PYTHONS['2.7mu'] = '/opt/python/cp27-cp27mu/bin/python'

    def __init__(self,
                 tarfile,
                 pytraj_home,
                 python_versions,
                 use_manylinux=False,
                 cpptraj_dir=''):
        self.libcpptraj = ''  # will be updated later
        self.is_osx = sys.platform.startswith('darwin')
        self.python_versions = python_versions
        self.use_manylinux = use_manylinux
        if cpptraj_dir:
            self.cpptraj_dir = os.path.abspath(cpptraj_dir)
        else:
            self.cpptraj_dir = ''
        self.miniconda_root = self.find_miniconda_root()
        if self.use_manylinux:
            self.python_exe_paths = dict(
                (py_version, PipBuilder.MANY_LINUX_PYTHONS[py_version])
                for py_version in self.python_versions)
        else:
            self.python_exe_paths = dict((py_version, '/'.join(
                (self.miniconda_root, '/envs/pytraj' + py_version,
                 'bin/python'))) for py_version in self.python_versions)
        self.repair_exe = ([
            pytraj_home + '/scripts/misc/fix_rpath_pip_osx.py',
        ] if self.is_osx else ['auditwheel', 'repair'])

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
        cmlist = '{python} -m pip wheel {tarfile}'.format(
            python=python, tarfile=tarfile).split()
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
        cmlist = 'conda create -n {} python={} numpy nomkl --yes'.format(
            env, python_version)
        print(cmlist)
        subprocess.check_call(cmlist.split())

    def _get_wheel_file(self, python_version, folder='./'):
        pdict = dict(darwin='macos', linux='linux')
        platform = pdict[sys.platform]
        whl_file = glob(folder + '/*' + python_version.replace('.', '') +
                        '*{}*.whl'.format(platform))[0]
        print(whl_file)
        return whl_file

    def repair_wheel(self, python_version):
        print('repair_wheel')
        whl_name = self._get_wheel_file(python_version)
        command = self.repair_exe + [
            whl_name,
        ]
        if self.is_osx:
            command.extend(['--libcpptraj', self.libcpptraj])
        print(command)
        subprocess.check_call(command)

    def _check_numpy_and_fix(self, python_exe, env):
        try:
            subprocess.check_call(
                '{} -c "import numpy"'.format(python_exe), shell=True)
        except subprocess.CalledProcessError:
            subprocess.check_call(
                '{} -m pip install numpy'.format(python_exe), shell=True)

    def validate_install(self, py_version):
        python_exe = self.python_exe_paths[py_version]
        env = 'pytraj' + py_version
        # e.g: change 2.7 to 27
        print('Testing pytraj build')
        whl_file = self._get_wheel_file(py_version, folder='wheelhouse')
        print('Testing wheel file {}'.format(whl_file))
        try:
            subprocess.check_call(
                '{} -m pip uninstall pytraj -y'.format(python_exe), shell=True)
        except subprocess.CalledProcessError:
            pass

        # always try with newest pip
        subprocess.check_call(
            '{} -m pip install pip --upgrade'.format(python_exe), shell=True)
        subprocess.check_call(
            '{} -m pip install {}'.format(python_exe, whl_file), shell=True)
        self._check_numpy_and_fix(python_exe, env)
        if sys.platform.startswith('darwin'):
            print('libcpptraj', self.libcpptraj)
            with temporarily_move_libcpptraj(self.libcpptraj):
                # moving libcpptraj to make sure pytraj use libcpptraj in pytraj/lib/
                # this is for osx only
                output = subprocess.check_output(
                    '{} -c "import pytraj as pt; pt.run_tests()"'.format(
                        python_exe),
                    shell=True)
        else:
            output = subprocess.check_output(
                '{} -c "import pytraj as pt; pt.run_tests()"'.format(
                    python_exe),
                shell=True)
        output = output.decode()
        print(output)
        print('PASSED: build test for {}'.format(whl_file))
        subprocess.check_call(
            '{} -m pip uninstall pytraj -y'.format(python_exe).split())

    def check_cpptraj_and_required_libs(self):
        pytraj_home = os.path.abspath(
            os.path.dirname(__file__).strip('scripts'))
        cpptraj_home = os.environ.get("CPPTRAJHOME", '')
        ext = 'so' if not sys.platform.startswith('darwin') else 'dylib'
        if not cpptraj_home:
            print('CPPTRAJHOME is not set')
            if not self.cpptraj_dir:
                self.cpptraj_dir = pytraj_home + '/cpptraj/'
            suggested_libcpptraj = self.cpptraj_dir + '/lib/libcpptraj.' + ext
            print('Looking for {}'.format(suggested_libcpptraj))
            if os.path.exists(suggested_libcpptraj):
                print('Found')
                os.environ['CPPTRAJHOME'] = self.cpptraj_dir
                print('CPPTRAJHOME is set to {}'.format(self.cpptraj_dir))
            else:
                print("Can not find libcpptraj")
                raise EnvironmentError(
                    'Must set CPPTRAJHOME having lib/libcpptraj')
        else:
            self.cpptraj_dir = cpptraj_home
        self.libcpptraj = os.path.join(self.cpptraj_dir, 'lib', 'libcpptraj.' + ext)

        for package in self.REQUIRED_PACKAGES:
            try:
                __import__(package)
            except ImportError:
                raise ImportError('require {}'.format(package))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Build wheel file to upload to pypi')
    parser.add_argument('tarfile')
    parser.add_argument(
        '--py',
        default=None,
        help='Python version. Default: build all versions (2.7, 3.4, 3.5)')
    parser.add_argument(
        '--cpptraj-dir', default='', help='cpptraj dir, optional')
    parser.add_argument(
        '--manylinux-docker',
        action='store_true',
        help='If specified, use Python versions from manylinux')
    args = parser.parse_args()

    tarfile = os.path.abspath(args.tarfile)
    python_versions = PipBuilder.SUPPORTED_VERSIONS if args.py is None else [
        args.py,
    ]
    # pytraj tar file
    pytraj_home = os.path.dirname(__file__).strip('scripts')
    builder = PipBuilder(
        tarfile=tarfile,
        pytraj_home=pytraj_home,
        python_versions=python_versions,
        use_manylinux=args.manylinux_docker,
        cpptraj_dir=args.cpptraj_dir)
    builder.run()
    # builder.libcpptraj = '../cpptraj/lib/libcpptraj.dylib'
    # builder.validate_install('3.5')
