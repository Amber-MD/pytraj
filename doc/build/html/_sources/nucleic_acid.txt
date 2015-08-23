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
    traj.top.residue_names
    na = pt.nastruct(traj)
    na
    na.mean_and_std()
    na.to_dict()
