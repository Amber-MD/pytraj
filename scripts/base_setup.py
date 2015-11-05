import os

message_auto_install = """
Can not find cpptraj header and libcpptraj files.
We're trying to dowload and build libcpptraj for you, would take about 5-10 minutes.
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
"""

message_openmp_cpptraj = '''
libcpptraj was detected to be installed with openmp.
You can not use --disable-openmp flag with pytraj
'''

message_serial_cpptraj = '''
libcpptraj was detected not be installed with openmp. You can recompile it with -openmp flag or 
disable openpm install in pytraj by adding --disable-openmp

Example: python setup.py install --disable-openmp
'''

message_after_sucessful_install = '''
make sure to add {0} to your LD_LIBRARY_PATH
    - Example: export LD_LIBRARY_PATH={1}:$LD_LIBRARY_PATH

Run test:
    - simple (few seconds): python ./runtests.py simple
    - full (5-10 minutes): python runtests.py

Check compiling info:
    # change directory so python does not try to import source code
    cd tests # or any other folders
    python -c 'import pytraj as pt; pt.show_versions()'

'''

def remind_export_LD_LIBRARY_PATH(build_tag, libdir, pytraj_inside_amber):
    if build_tag:
        if not pytraj_inside_amber:
            from scripts.acsii_art import batman
            print("")
            print("")
            print(batman)
            libdir = os.path.abspath(libdir)
            print(message_after_sucessful_install.format(libdir, libdir))
            print("")
        else:
            # pytraj is a part of Amber
            print('make sure to `source $AMBERHOME/amber.sh` (if using bash) '
                  'or `source $AMBERHOME/amber.csh` if using csh')
    else:
        print("not able to install pytraj")
