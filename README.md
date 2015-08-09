[![Build Status](https://travis-ci.org/Amber-MD/pytraj.svg?branch=master)](https://travis-ci.org/Amber-MD/pytraj)
[![Binstar Badge](https://binstar.org/pytraj/pytraj-dev/badges/version.svg)](https://binstar.org/pytraj/pytraj-dev/)

PYTRAJ
------

- pytraj is a Python front-end to cpptraj program (a data analysis for biomolecular simulation)
- Why using pytraj?
    * It's fast
        * most of `pytraj's codes` were written in [Cython language] (http://cython.org/)
        * its core (cpptraj) was written in C++ (big thanks to [cpptraj developers] (https://github.com/mojyt/cpptraj))
        * it supports parallel computing (openmp from cpptraj or mpi from mpi4py or parallel in ipython)
    * It has both in-memory and out-of-core calculations
    * It supports more than 100 types of analyses in [cpptraj] (http://ambermd.org/doc12/Amber15.pdf)
    * It has clean Python/Cython syntax
    * It has been extensively tested (>10K lines of testing code)
    * It's portable: you only need to install "libcpptraj" and optional numpy
- Note: `pytraj` is still in its infancy and its API might be rapidly changed. But it's not hurt to try :).


Install
-------
- require:
    - cpptraj
    - numpy
- (optional):
    - libnetcdf (highly recommended)
    - pandas (if you want to use DataFrame) 
    - matplotlib (if you want to plot)
- easiest and less headache:
    * git clone https://github.com/pytraj/pytraj
    * cd pytraj
    * python ./setup.py install
    * (Note: pytraj will install the most update cpptraj)
- if you are using `conda`, you can just `conda install -c pytraj pytraj-dev` for Linux system
- further: check wiki page [wiki](http://www.github.com/pytraj/pytraj/wiki)

How to get started?
------------------
- examples: 

    ```python
    >>> import pytraj as pt
    >>> traj = pt.iterload("data.nc", "top.parm7")
    >>> pt.rmsd(traj, ref=0, mask='@CA')
    >>> pt.distance(traj, [[0, 2], [3, 7]])
    >>> pt.bfactors(traj, '@CA', byres=True, dtype='dataset').plot()
    >>> pt.energy_decomposition(traj, igb=8, parm="./top.parm7")['dihedral']
    >>> pt.dssp(traj, ':2-16')
    >>> pt.calc_phi(traj, resrange=range(2, 8, 2))
    >>> pt.rotate_dihedral(traj, '3:phi:120')
    >>> traj['@CA'].xyz

    ```
- many more:
    * check [examples folder] (./examples/) (always keep updated)
    * check pytraj-notebook: [pytraj-notebook](http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/pytraj_overview.ipynb) (sometimes outdated)
    * [more will come] (http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/index.ipynb)
    * pytraj's tutorials and documents are growing, stay tuned.

Question/Suggestion?
--------------------
* code issue and stuff relating directly to `pytraj`, create [Issue](https://github.com/pytraj/pytraj/issues)
* ask question about data analysis in general, send email to [AMBER Mailing List] (http://lists.ambermd.org/mailman/listinfo/amber)

Support
-------
* Development version of [cpptraj] (https://github.com/mojyt/cpptraj)

[Contributors and Acknowledgement] (./contributors/)
----------------------------------------------------

Citation
--------
- cpptraj : [PTRAJ and CPPTRAJ] (http://pubs.acs.org/doi/abs/10.1021/ct400341p): Software for Processing and Analysis of Molecular Dynamics Trajectory Data
Daniel R. Roe and Thomas E. Cheatham, III
Journal of Chemical Theory and Computation 2013 9 (7), 3084-3095 

- pytraj : Hai Nguyen et al. (2015) (in preperation)

License
-------
BSD 2-Clause
