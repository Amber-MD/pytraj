[![Build Status](https://travis-ci.org/hainm/pytraj.svg?branch=master)](https://travis-ci.org/hainm/pytraj)
Welcome to pytraj!
-------------------

- pytraj is a Python package wrapping cpptraj program (a data analysis for biomolecular simulation)
- Why using pytraj:
    * It's fast
        * its core (cpptraj) was written in C++
        * it has interface with numpy. Data calculation are performed without copying to numpy array
        * (but it actually does not need `numpy` at all, just optional)
    * It has clean Python/Cython syntax
    * It has been extensively tests
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
    - netcdf library
    - cpptraj
- wiki page : [wiki](http://www.github.com/pytraj/pytraj/wiki)
- install libcpptraj: 
    ./installs/libcpptraj.txt (works well with development version of cpptraj)
- installs pytraj: [wiki](http://www.github.com/pytraj/pytraj/wiki)
    (further instruction ./installs/pytraj.txt)

Usage: 
-----
- easiest: dist = `calculate('distance', ':2@CA :10@CA', (traj, traj))`
-                 `calculate('matrix', '@CA', frame, top)`
- Check ./examples folder
- Check pytraj-notebooks: [pytraj-notebook](http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/Frame_class.ipynb)

Support
====================
* Development version of [cpptraj] (https://github.com/mojyt/cpptraj)
