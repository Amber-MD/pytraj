[![Build Status](https://travis-ci.org/Amber-MD/pytraj.svg?branch=master)](https://travis-ci.org/Amber-MD/pytraj)
[![Binstar Badge](https://binstar.org/ambermd/pytraj-dev/badges/version.svg)](https://binstar.org/ambermd/pytraj-dev/)
[![Coverage Status](https://coveralls.io/repos/Amber-MD/pytraj/badge.svg?branch=master&service=github)](https://coveralls.io/github/Amber-MD/pytraj?branch=master)

Try online [![Binder](http://mybinder.org/images/logo.svg)](http://mybinder.org/repo/hainm/notebook-pytraj)
-----------------------------------------------------------------------------------------------------------

PYTRAJ
------

pytraj is a Python front-end of cpptraj program (a data analysis package for biomolecular simulation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

http://amber-md.github.io/pytraj

[![pytraj website](./examples/progress_bar.gif)](http://amber-md.github.io/pytraj/latest/index.html)


Install
-------

- from source:
    * git clone https://github.com/amber-md/pytraj
    * cd pytraj
    * python ./setup.py install
    * (Note: pytraj will install the most updated cpptraj)
- from conda (linux onlye): `conda install -c ambermd pytraj-dev libcpptraj-dev`
- getting trouble? : [check our webpage](http://amber-md.github.io/pytraj/doc/build/html/installation.html)

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

nglview with pytraj in Jupyter notebook
---------------------------------------

[![pytraj website](./examples/figures/nglview_pytraj.gif)](http://amber-md.github.io/pytraj/latest/index.html)

License
-------
BSD 2-Clause

