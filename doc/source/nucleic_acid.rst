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

    # get values for major groove width
    na.major

    # get values for several parameters
    for key in ['minor', 'twist', 'incl']:
        print(na[key])
