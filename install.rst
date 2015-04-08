# Install
## Using AmberTools 15 (available after April 2015): [AmberTools 15] (http://ambermd.org/#AmberTools)
    * if you are familiar with AMBER and AmberTools, you can install `pytraj` from:
          * cd $AMBERHOME
          * ./configure gnu
          * cd $AMBERHOME/AmberTools/src/
          * make pytraj
          * (for further info, check `Installation` section in Amber manual
          * Question? Send email to [AMBER mailing-list]  (http://ambermd.org/#reflector)

## Using[conda] (http://conda.pydata.org/docs/examples/install.html)
    * conda install -c pytraj pytraj-dev --force

## Using github version
* if using the most updated pytraj version (recommended)
    * `git clone https://github.com/pytraj/pytraj`
    * `cd pytraj`
    * `python ./setup.py install`
* if using specific verion (go to [releases](https://github.com/pytraj/pytraj/releases/) to check the most updated and stable version)
    * `wget https://github.com/pytraj/pytraj/archive/v0.1.2.dev0.1.tar.gz`
    * `tar -xf v0.1.2.dev0.1.tar.gz`
    * `cd pytraj*`
    * `python setup.py install`

## Using pip
* Note: There is no development version on pypi.
* install libcpptraj
    * `git clone https://github.com/mojyt/cpptraj`
    * `cd cpptraj`
    * export CPPTRAJHOME=\`pwd\`
    * `./configure -shared -nomathlib gnu`
        * (if you're using Amber, you can do `./configure -shared -nomathlib -amberlib gnu`)
    * `make libcpptraj`
* then install pytraj
    * (make sure to `export CPPTRAJHOME=your_cpptraj_dir` because `pytraj` will find header files in `$CPPTRAJHOME/src` and libcpptraj in `$CPPTRAJHOME/lib`
    * `pip install pytraj`

# Ruu small test after installing
`python -c 'import pytraj as pt; pt.run_tests()'`

if you install via github, you will get full test runs
* `cd tests`
* `python ./run_all_and_find_fails.py

Possible errors:
> 1. "ImportError: libcpptraj.so: cannot open shared object file: No such file or directory"

> > you need to add libcpptraj.so to LD_LIBRARY_PATH
