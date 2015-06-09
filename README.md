[![Build Status](https://travis-ci.org/pytraj/pytraj.svg?branch=master)](https://travis-ci.org/pytraj/pytraj)
[![Binstar Badge](https://binstar.org/pytraj/pytraj-dev/badges/version.svg)](https://binstar.org/pytraj/pytraj-dev/)

PYTRAJ
------

- pytraj is a Python package wrapping cpptraj program (a data analysis for biomolecular simulation)
- Why using pytraj?
    * It's fast
        * most of `pytraj's codes` were written in [Cython language] (http://cython.org/)
        * its core (cpptraj) was written in C++ (big thanks to [cpptraj developers] (https://github.com/mojyt/cpptraj))
        * it supports parallel computing (openmp from cpptraj or mpi from mpi4py or parallel in ipython)
        * it has interface with numpy. Data calculation are performed without copying
    * It supports more than 100 types of analyses in [cpptraj] (http://ambermd.org/doc12/Amber15.pdf)
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


Install
-------
- require:
    - cpptraj
- (optional):
    - libnetcdf (highly recommended)
    - numpy (highly recommended)
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

How to get started?
------------------
- examples: 
    * calculate distance: `dist = calc_distance(traj, ':2@CA :10@CA')`
    * calculate distance matrix: `mat = calc_matrix(frame, '@CA', top)`
    * calculate DSSP: `calc_dssp([[frame,], traj1, traj2(3, 9, 2), traj3.chunk_iter(chunk=5)], ':2-10', dtype='string', top=traj.top)`
    * get new Trajectory with a given mask: `traj['@CA']`
    * update coords for specific atoms in residues: `traj[':1-3'] = xyz_3d`
    * scale and translate coords for non-H atoms: `traj.apply(lambda x : x * 2 + 1.2, indices_or_mask='!@H='`)
    * expose to numpy: `arr0 = np.asarray(frame)` 
    * convert to pandas's DataFrame: `dframe = traj.search_hbonds(dtype='dataframe')`
    * load from other packages: 
        * from `mdtraj`: `traj = io.load_mdtraj(mdtraj_traj_object)`
        * from `MDAnalysis`: `trajiterator = io.load_MDAnalysisIterator(universe_object)`
        * from `ParmEd`: `parm = io.load_ParmEd(its_obj)`
    * expose to Cython (will be translated to C++ code): `from pytraj.Frame cimport _Frame`
- many more:
    * check ./examples folder
    * check pytraj-notebook: [pytraj-notebook](http://nbviewer.ipython.org/github/pytraj/pytraj/blob/master/note-books/pytraj_overview.ipynb)
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
