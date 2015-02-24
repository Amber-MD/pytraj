# Install
## 1. github version (pytraj will take care of installing `libcpptraj`)
* if using the most updated pytraj version (recommended)
    * `git clone https://github.com/pytraj/pytraj`
    * `cd pytraj`
    * `python ./setup.py install`
* if using specific verion (go to [releases](https://github.com/pytraj/pytraj/releases/) to check the most updated and stable version)
    * `wget https://github.com/pytraj/pytraj/archive/0.1.beta.8.tar.gz`
    * `tar -xf 0.1.beta.8.tar.gz`
    * `cd pytraj-0.1.beta.8`
    * `python setup.py install`
* Done
* Simple test: python -c 'import pytraj as pt; pt.run_tests()'

> Note: Possible failure: can not find netcdf.h file. Make sure to install netcdeflib or specify its include path

## 2. Using pip
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

## 3. Using conda
Not yet

## 4. run small test after installing
`python -c 'import pytraj as pt; pt.run_tests()'`

if you install via github, you will get full test runs
* `cd tests`
* `python ./RunAllAndFindFailure.py`

Possible errors:
> 1. "ImportError: libcpptraj.so: cannot open shared object file: No such file or directory"

> > you need to add libcpptraj.so to LD_LIBRARY_PATH
