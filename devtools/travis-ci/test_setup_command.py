#!/usr/bin/env python
import os
import subprocess
from nose import tools

# tested with py3.5

def capture_stderr(command):
    try:
        subprocess.check_output(command.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        return e.output.decode()

def test_default_raise_if_install_cpptraj_without_openmp():
    """ install libcpptraj without openmp """
    subprocess.check_output('python ./scripts/install_libcpptraj.py'.split())
    command = 'python setup.py install'
    output = capture_stderr(command)
    tools.assert_in('libcpptraj was NOT installed with openmp', output)
    tools.assert_in('Turn off openmp in pytraj: python setup.py install --disable-openmp', output)

def test_setup_clean():
    """ python setup.py clean """
    command = 'python setup.py clean'
    output = subprocess.check_output(command.split()).decode()
    tools.assert_in('must_compile_c_extension =  False', output)
    tools.assert_in('use_pip = False', output)
    tools.assert_in('running clean', output)
    tools.assert_not_in('libcpptraj', output)
    tools.assert_not_in('github', output)

def test_setup_help():
    """ python setup.py --help """
    for command in ['python setup.py --help', 'python setup.py -h']:
        output = subprocess.check_output(command.split()).decode()
        tools.assert_in('setup.py install    will install the package', output)
        tools.assert_not_in('libcpptraj', output)
        tools.assert_not_in('github', output)

def test_pip_help():
    """ pip install --help """
    command = 'pip install --help'
    output = subprocess.check_output(command.split()).decode()
    tools.assert_in('Install Options', output)

def test_raise_if_using_pip_but_does_not_have_cpptraj_home():
    os.environ['CPPTRAJHOME'] = ''
    command = 'pip install -e .'
    output = capture_stderr(command)
    tools.assert_in('must_compile_c_extension =  True', output)
    print(output)
