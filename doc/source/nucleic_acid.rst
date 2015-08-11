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
    traj.top.atomnames
    na = pt.nastruct(traj)
    na
    na.summary()
    na.values
