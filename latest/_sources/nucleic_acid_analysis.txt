.. _nucleic_acid:

Nucleic Acid Analysis
=====================

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)

.. autofunction:: pytraj.nastruct 

Longer Examples
---------------

.. ipython:: python
 
    import pytraj as pt
    traj = pt.iterload('data/adh026.3.pdb')
    traj

    # check unique residue names
    print(set(residue.name for residue in traj.top.residues))
    na = pt.nastruct(traj)
    na

    # get values for major groove width
    na.major

    # get values for several parameters
    for key in ['minor', 'twist', 'incl']:
        print(na[key])

    # get some statistics (syntax might be changed in future)
    import numpy as np
    # get mean
    na._summary(np.mean)

    # get std
    na._summary(np.std)

    # only interested in some parameters?
    na._summary(np.mean, keys=['major', 'minor', 'twist'])

    # explain keywords
    print(na._explain())

    # if we have long analysis, we can temporarily save ``na`` to disk by `pytraj.to_pickle`` and load back later.
    pt.to_pickle(na, 'na_.pk')

    # load back pickle object
    na2 = pt.read_pickle('na_.pk')
    na2
    na2._summary(np.mean, ['major', 'twist'])
