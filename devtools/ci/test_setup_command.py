#!/usr/bin/env python
import os
from contextlib import contextmanager
import subprocess


# tested with py3.6
# utils
def clone_cpptraj():
    os.system('git clone https://github.com/amber-md/cpptraj')


@contextmanager
def temporarily_rename_cpptraj_folder():
    old_name = 'cpptraj'
    new_name = 't_cpptraj'
    if os.path.exists(old_name):
        print(os.path.abspath(old_name))
        subprocess.check_call('mv {} {}'.format(old_name, new_name).split())
        yield
        subprocess.check_call('mv {} {}'.format(new_name, old_name).split())
    else:
        print("does not have cpptraj, nothing done")
        yield
        print("nothing done")


def capture_stderr(command):
    """capture_stderr"""
    try:
        subprocess.check_output(command.split(), stderr=subprocess.STDOUT)
        return None
    except subprocess.CalledProcessError as e:
        return e.output.decode()


def git_clean_folder(name):
    """git_clean_folder"""
    command = 'git clean -fdx {}'.format(name)
    print('clean: ', command)
    subprocess.check_call(command.split())


def get_output(command):
    return subprocess.check_output(command.split()).decode()


# do not need cpptraj folder
def test_setup_clean():
    """ python setup.py clean """
    command = 'python setup.py clean'
    output = subprocess.check_output(command.split()).decode()
    assert 'compile_c_extension = False' in output
    assert 'use_pip = False' in output
    assert 'running clean' in output
    assert 'libcpptraj' not in output
    assert 'github' not in output


def test_setup_help():
    """ python setup.py --help """
    for command in ['python setup.py --help', 'python setup.py -h']:
        output = subprocess.check_output(command.split()).decode()
        assert 'setup.py install    will install the package' in output
        assert 'libcpptraj' not in output
        assert 'github' not in output


def test_pip_help():
    """ pip install --help """
    command = 'pip install --help'
    output = subprocess.check_output(command.split()).decode()
    assert 'Install Options' in output


def test_raise_if_using_pip_but_does_not_have_cpptraj_home():
    """ test_raise_if_using_pip_but_does_not_have_cpptraj_home_or_not_cpptraj_in_this_folder
        pip install -e .
        pip install .
        pip wheel .
    """
    os.environ['CPPTRAJHOME'] = ''
    with temporarily_rename_cpptraj_folder():
        for command in ['pip install -e .', 'pip install .', 'pip wheel .']:
            print(command)
            output = capture_stderr(command)
            assert 'Detected use of pip; you must set CPPTRAJHOME' in output


def test_install_libcpptraj_if_having_cpptraj_folder_here():
    clone_cpptraj()
    git_clean_folder('./cpptraj/lib/libcpptraj*')
    command = 'python setup.py build'
    output = get_output(command)
    assert 'cpptraj/lib/"libcpptraj' in output
    assert 'install libcpptraj from current' in output


def test_install_to_amberhome():
    fn = './fake_amberhome'
    python_path = '{}/lib/python3.6/site-packages/'.format(fn)
    mkdir_cmd = 'mkdir -p {}'.format(python_path)
    subprocess.check_call(mkdir_cmd, shell=True)
    full_name = os.path.abspath(fn)
    command = 'PYTHONPATH={} python setup.py install --prefix={}'.format(
        python_path,
        full_name)
    print(command)
    output = subprocess.check_output(command, shell=True).decode()
    assert 'Finished processing dependencies for pytraj' in output
    git_clean_folder(full_name)


def test_default_raise_if_install_cpptraj_without_openmp():
    """ test_default_raise_if_install_cpptraj_without_openmp """
    subprocess.check_output('python ./scripts/install_libcpptraj.py'.split())
    command = 'python setup.py install'
    output = capture_stderr(command)
    assert 'libcpptraj was NOT installed with openmp' in output
    assert 'Turn off openmp in pytraj: python setup.py install --disable-openmp' in output
