[![Build Status](https://travis-ci.org/Amber-MD/pytraj.svg?branch=master)](https://travis-ci.org/Amber-MD/pytraj)
[![Binstar Badge](https://binstar.org/ambermd/pytraj-dev/badges/version.svg)](https://binstar.org/ambermd/pytraj-dev/)
[![Coverage Status](https://coveralls.io/repos/Amber-MD/pytraj/badge.svg?branch=master&service=github)](https://coveralls.io/github/Amber-MD/pytraj?branch=master)
[![Binder](http://mybinder.org/images/logo.svg)](http://mybinder.org/repo/hainm/notebook-pytraj)

PYTRAJ
------
A Python front-end of [cpptraj program] (https://github.com/Amber-MD/cpptraj) (a data analysis package for biomolecular simulation): http://amber-md.github.io/pytraj

* [Features](#features)
* [Install](#install)
* [How to get started?](#how-to-get-started)
* [Contributors and Acknowledgement](#contributors-and-acknowledgement)
* [Citation](#citation)
* [Question/Suggestion?](#questionsuggestion)
* [nglview with pytraj in Jupyter notebook](#nglview-with-pytraj-in-jupyter-notebook)
* [Demo: Interactive data exploration with Jupyter notebook](#demo-interactive-data-exploration-with-jupyter-notebook)
* [License](#license)


Features
--------

- support more than 80 types of data analyses (rmsd, radgyr, autoimage, pca, clustering,...)
- read/write various file formats (.nc, .dcd, .trr, .pdb, .mol2, ...)
- fast (core codes were written in C++ and Cython)
- support parallel calculation with trivial installation (openmp, multiprocessing, mpi, ...)
- interactive analysis with large trajectory data that does not fit to memory
- [>> many more with comprehensive tutorials] (http://amber-md.github.io/pytraj)


Install
-------

- from source:

    ```bash
    git clone https://github.com/amber-md/pytraj
    cd pytraj

    # linux
    python ./setup.py install
 
    # osx: not well tested (yet)
    python ./setup.py install --disable-openmp
    ```

- from conda (develop version for linux, osx): `conda install -c https://conda.anaconda.org/ambermd pytraj-dev`

- from AMBER suite distribution [http://ambermd.org/](http://ambermd.org/). AMBER16 will be available soon (around April 2016).

- getting trouble? : [check our webpage](http://amber-md.github.io/pytraj/latest/installation.html)


How to get started?
------------------

- examples: 

    ```python
    import pytraj as pt
    traj = pt.iterload("data.nc", "top.parm7")
    pt.rmsd(traj, ref=0, mask='@CA')
    pt.dssp(traj, ':2-16')
    ```
- check our website: [http://amber-md.github.io/pytraj] (http://amber-md.github.io/pytraj)

Contributors and Acknowledgement
--------------------------------

[Please check here](./contributors/)

Citation
--------

If you would like to acknowledge our work, please cite both ``cpptraj`` and ``pytraj`` papers

- [PTRAJ and CPPTRAJ] (http://pubs.acs.org/doi/abs/10.1021/ct400341p): Software for Processing and Analysis of Molecular Dynamics Trajectory Data
Daniel R. Roe and Thomas E. Cheatham, III
Journal of Chemical Theory and Computation 2013 9 (7), 3084-3095 

- PYTRAJ: Interactive data analysis for molecular dynamics simulations. Hai Nguyen, Daniel R. Roe, Jason Swails, David A. Case. (2016)

Question/Suggestion?
--------------------
* code issue and stuff relating directly to `pytraj`, create [Issue](https://github.com/pytraj/pytraj/issues)
* ask question about data analysis in general, send email to [AMBER Mailing List] (http://lists.ambermd.org/mailman/listinfo/amber)

[nglview](https://github.com/arose/nglview) with pytraj in Jupyter notebook
---------------------------------------

[![pytraj website](./examples/figures/nglview.png)](http://amber-md.github.io/pytraj/latest/index.html)

Demo: Interactive data exploration with [Jupyter notebook] (http://jupyter.org/)
--------------------------------------------------------------------------------

[![pytraj website](./examples/progress_bar.gif)](http://amber-md.github.io/pytraj/latest/index.html)


License
-------
BSD 2-Clause
