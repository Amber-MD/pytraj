.. _nucleic_acid:

Nucleic Acid Analysis
=====================

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)

.. autofunction:: pytraj.nastruct 

Examples
--------

.. ipython:: python
 
    import pytraj as pt
    traj = pt.iterload('data/adh026.3.pdb')
    traj

    # check unique residue names
    traj.top.residue_names
    na = pt.nastruct(traj)
    na

    # get mean and std for each analysis (major groove width, minor, ...)
    na.mean_and_std()

    # get raw data for major groove
    na.dataset.filter(lambda x : 'major' in x.key)
