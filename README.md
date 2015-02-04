[![Build Status](https://travis-ci.org/hainm/pytraj.svg?branch=master)](https://travis-ci.org/hainm/pytraj)
Welcom to pytraj!
-------------------

- pytraj is a Python package wrapping cpptraj program (a data analysis for biomolecular simulation)
- Why using pytraj:
    * It's fast
        * its core (cpptraj) was written in C++
        * it has interface with numpy. Data calculation are performed without copying to numpy array
        * (but it actually does not need `numpy` at all, just optional)
    * It has clean Python/Cython syntax
    * It has been extensively testes
    * It's flexible: 
        * you can write extension modules in either high (Python) or low (C/C++ or Cython) level

Citing (optional):
-----------------
- cpptraj : [PTRAJ and CPPTRAJ] (http://pubs.acs.org/doi/abs/10.1021/ct400341p): Software for Processing and Analysis of Molecular Dynamics Trajectory Data
Daniel R. Roe and Thomas E. Cheatham, III
Journal of Chemical Theory and Computation 2013 9 (7), 3084-3095 

- pytraj : (in preparation)

Install
-------
- install libcpptraj: 
    ./installs/libcpptraj.txt (works well with cpptraj v15.22b)
- installs pytraj: `pip install pytraj` or `pip install --upgrade pytraj`
    (further instruction ./installs/pytraj.txt)
- wiki page : [wiki](http://www.github.com/pytraj/pytraj/wiki)

Usage: 
-----
- Check ./examples/ folders
- Chec pytraj-notebooks: [pytraj-notebook](http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/Frame_class.ipynb)

Version I am working on:
====================
* Development version of `cpptraj`: v15.22b
* Python 2.7.8 :: Anaconda 2.1.0 (64-bit)
* Cython 0.22pre (pre-verion of Cython 0.22)
    * Development version: https://github.com/cython/cython
