[![Build Status](https://travis-ci.org/hainm/pytraj.svg?branch=master)](https://travis-ci.org/hainm/pytraj)
[![Binstar Badge](https://binstar.org/pytraj/pytraj-dev/badges/version.svg)](https://binstar.org/pytraj/pytraj-dev/)
Welcome to pytraj!
-------------------
[wiki](http://www.github.com/pytraj/pytraj/wiki)

- pytraj is a Python package wrapping cpptraj program (a data analysis for biomolecular simulation)
- Why using pytraj:
    * It's fast
        * its core (cpptraj) was written in C++
        * it has interface with numpy. Data calculation are performed without copying to numpy array
        * (but it actually does not need `numpy` at all, just optional)
    * It has clean Python/Cython syntax
    * It has been extensively tested 
    * It's flexible: 
        * you can write extension modules in either high (Python) or low (C/C++ or Cython) level

Citation:
-----------------
- cpptraj : [PTRAJ and CPPTRAJ] (http://pubs.acs.org/doi/abs/10.1021/ct400341p): Software for Processing and Analysis of Molecular Dynamics Trajectory Data
Daniel R. Roe and Thomas E. Cheatham, III
Journal of Chemical Theory and Computation 2013 9 (7), 3084-3095 

- pytraj : Hai Nguyen et al. (2015) (in preperation)

Install
-------
- require:
    - cpptraj
    - netcdf library (optional but highly recommended)
- easiest:
    * git clone https://github.com/pytraj/pytraj
    * cd pytraj
    * python ./setup.py install
    * (It's really slow to install? Try `CFLAGS="-O0 -ggdb" python ./setup.py install`)
- further: check wiki page : [wiki](http://www.github.com/pytraj/pytraj/wiki)
- install libcpptraj: 
    ./installs/libcpptraj.txt (works well with development version of cpptraj)
- installs pytraj: [wiki](http://www.github.com/pytraj/pytraj/wiki)
    (further instruction ./installs/pytraj.txt)
- if you are using `conda`, you can just `conda install -c pytraj pytraj-dev` for Linux system

Usage: 
-----
- example: 
    * `dist = calculate('distance', ':2@CA :10@CA', (traj, traj))`
    * `mat = calculate('matrix', '@CA', frame, top)`
    * get 2D array (xyz coords)for a given mask: `traj['@CA']`
    * expose to numpy: `arr0 = np.asarray(frame[:])` 
    * load from other package: `traj = FrameArray(mdtraj_traj.xyz, top)`
    * expose to Cython (will be translated to C++ code): `from pytraj.Frame cimport _Frame`
- many more:
    * check ./examples folder
    * check pytraj-notebook: [pytraj-notebook](http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/pytraj_overview.ipynb)
    * [more will come] (http://nbviewer.ipython.org/github/pytraj/pytraj/tree/master/note-books/)
    * pytraj's document is still poorly written, so you can also check ./tests folder for more 

Question?
--------
* send email to [AMBER Mailing List] (http://lists.ambermd.org/mailman/listinfo/amber)

Support
====================
* Development version of [cpptraj] (https://github.com/mojyt/cpptraj)
