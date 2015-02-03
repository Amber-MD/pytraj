[![Build Status](https://travis-ci.org/hainm/pytraj.svg?branch=master)](https://travis-ci.org/hainm/pytraj)
Welcom to pytraj!
-------------------


Aim: Write cpptraj wrapper for Python user
-----------------------------------------
- Short term goal: expose cpptraj library to Python
- Very long term goal: glue most of useful library in [AMBER] (http://ambermd.org/) by Python: 
    * tleap --> pyleap
    * paramfit --> pyparamfit
    * reduce --> pyreduce 
    * ParmedTools
    * (add many more)
- Why using pytraj:
    * It's fast
        * it's a wrapper of cpptraj (was written in C++ (by Daniel R. Roe))
        * it has interface with numpy. Data calculation are performed without copying to numpy array
        * (but it actually does not need `numpy` at all, just optional)
    * It has clean syntax
        * Python/Cython style with fancy indexing 
    * It has been extensively testes
    * It's flexible: 
        * you can write extension modules in either high (Python) or low (C/C++ or Cython) level

Status: pre-release version 0.1
------------------------------

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
- wiki page : [wiki](http://www.github.com/hainm/pytraj/wiki)

Usage: 
-----
- Make sure to export LD_LIBRARY_PATH before using
    + export LD_LIBRARY_PATH=$AMBERHOME/lib/:$LD_LIBRARY_PATH
    + or export LD_LIBRARY_PATH=$CPPTRAJHOME/lib:$LD_LIBRARY_PATH

- Check ./examples/

Version I am working on:
====================
* Development version of `cpptraj`: v15.22b
* Python 2.7.8 :: Anaconda 2.1.0 (64-bit)
* Cython 0.22pre (pre-verion of Cython 0.22)
    * Development version: https://github.com/cython/cython
