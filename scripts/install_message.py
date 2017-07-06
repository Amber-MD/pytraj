__all__ = [
    'message_cython', 'message_pip_need_cpptraj_home', 'message_auto_install',
    'message_openmp_cpptraj', 'message_serial_cpptraj',
    'message_after_sucessful_install'
]

message_cython = '''
Warning:  pytraj was not installed !

Building from source requires cython >= 0.21

Either try:
    conda install cython
    pip install cython --upgrade

We suggest to install Anaconda suite, which has more than 300 python packages (including
the most updated cython)

    http://conda.pydata.org/docs/download.html)

'''

message_pip_need_cpptraj_home = '''

installing from pip for python 2.7 in OSX
requires to pre-install libcpptraj and to set CPPTRAJHOME

An example of installing libcpptraj:

$ git clone https://github.com/Amber-MD/cpptraj/
$ cd cpptraj
$ export CPPTRAJHOME=`pwd`
$ ./configure -shared -openmp gnu
$ make libcpptraj -j8

How to install pytraj within seconds without installing libcpptraj by yourself?

Either:
    - Use python >= 3.4 (OSX): pip install pytraj
    - Use conda (Linux, osx): conda install -c ambermd pytraj
'''

message_auto_install = '''
Can not find cpptraj header and libcpptraj files.
We're trying to download and build libcpptraj for you, would take about 5-10 minutes.
You can check ./cpptraj/ folder after installation.

To avoid auto-installation
--------------------------
You should set CPPTRAJHOME or installing ./cpptraj/ in this current folder.

If you want to manually install `libcpptraj`, you can download cpptraj
development version from here: https://github.com/Amber-MD/cpptraj

$ git clone https://github.com/Amber-MD/cpptraj/
$ cd cpptraj
$ export CPPTRAJHOME=`pwd`
$ ./configure -shared -openmp gnu
$ make libcpptraj

and then go back to pytraj folder:
python setup.py install
'''

message_openmp_cpptraj = '''
libcpptraj was detected to be installed with openmp.
You can not use --disable-openmp flag with pytraj
'''

message_serial_cpptraj = '''
libcpptraj was NOT installed with openmp. You can recompile it with -openmp flag or
disable openmp install in pytraj by adding --disable-openmp

Note: For osx users, pytraj uses clang for compiling cython extension and '--disable-openmp' flag
must be specified. If experienced users want to hack, please check setup.py file.

Examples:
    - Turn off openmp in pytraj: python setup.py install --disable-openmp
    - Turn ON openmp in cpptraj: ./configure -shared -openmp gnu && make libcpptraj -j4
    - Turn ON openmp in cpptraj and using netcdf, lapack, ... in $AMBERHOME (if having
      one):
          ./configure -shared -openmp -amberlib gnu && make libcpptraj -j4
'''

message_after_sucessful_install = '''
make sure to add {0} to your LD_LIBRARY_PATH
    - Example: export LD_LIBRARY_PATH={1}:$LD_LIBRARY_PATH

    - Notes: you can move `libcpptraj.so` to any where you want, just properly add it to $LD_LIBRARY_PATH

Run test:
    - simple (few seconds): python ./run_tests.py simple
    - full (5-10 minutes): python run_tests.py
'''
