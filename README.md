[![Build Status](https://travis-ci.org/hainm/pytraj.svg?branch=master)](https://travis-ci.org/hainm/pytraj)
[![Binstar Badge](https://binstar.org/pytraj/pytraj-dev/badges/version.svg)](https://binstar.org/pytraj/pytraj-dev/)
Welcome to pytraj!
-------------------
[wiki](http://www.github.com/pytraj/pytraj/wiki)

- pytraj is a Python package wrapping cpptraj program (a data analysis for biomolecular simulation)
- Why using pytraj?
    * It's fast
        * its core (cpptraj) was written in C++ with more than 80K lines of code (big thanks to [cpptraj developers] (https://github.com/mojyt/cpptraj))
        * it supports parallel computing (openmp from cpptraj or mpi from mpi4py or parallel in ipython)
        * it has interface with numpy. Data calculation are performed without copying to numpy array
        * (but it actually does not need `numpy`)
    * It supports more than 100 kinds of actions/analyses in [cpptraj] (http://ambermd.org/doc12/Amber14.pdf)
    * It has clean Python/Cython syntax
    * It has been extensively tested (>10K lines of testing code)
    * It's flexible: 
        * user does not need to care about file ext, `cpptraj` (and then `pytraj`) automatically detects it, as long as the file is supported
            (io.load("myparm.crazy_ext"))
        * you can write extension modules in either high (Python) or low (C/C++ or Cython) level
        * you can easily load objects from other packages (ParmEd, MDAnalysis, mdtraj...)
    * It's portable: you only need to install "libcpptraj"
        * (but you can use extra help from other popular packages such as numpy, matplotlib)
- Note: `pytraj` is still in its infancy and its API might be rapidly changed. But it's not hurt to try :).

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
- (optional):
    - libnetcdf (highly recommended)
    - numpy (just for interfacing with other packages since they use numpy)
    - pandas (if you want to use DataFrame) 
    - matplotlib (if you want to plot)
- easiest and less headache:
    * git clone https://github.com/pytraj/pytraj
    * cd pytraj
    * python ./setup.py install
        * if it's really slow to install? Try building in parallel
            * python ./setup.py build -faster_build
            * python ./setup.py install)
- further: check wiki page [wiki](http://www.github.com/pytraj/pytraj/wiki)
- if you are using `conda`, you can just `conda install -c pytraj pytraj-dev` for Linux system

Usage: 
-----
- example: 
    * `dist = calc_distance(traj, ':2@CA :10@CA')`
    * `mat = calc_matrix(frame, '@CA', top)`
    * `calc_dssp([[frame,], traj1, traj2(3, 9, 2), traj3.chunk_iter(chunk=5)], ':2-10', dtype='ndarray')`
    * get new Trajectory with a given mask: `traj['@CA']`
    * expose to numpy: `arr0 = np.asarray(frame[:])` 
    * load from other package: `traj = io.load_mdtraj(mdtraj_traj)`, `parm = io.load_ParmEd(its_obj)`
    * expose to Cython (will be translated to C++ code): `from pytraj.Frame cimport _Frame`
- many more:
    * check ./examples folder
    * check pytraj-notebook: [pytraj-notebook](http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/pytraj_overview.ipynb)
    * [more will come] (http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/index.ipynb)
    * pytraj's tutorials and documents are growing, stay tuned.

Question?
--------
* send email to [AMBER Mailing List] (http://lists.ambermd.org/mailman/listinfo/amber)

Support
====================
* Development version of [cpptraj] (https://github.com/mojyt/cpptraj)
